### training ###
import os
import sys
import copy

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

import comet_ml
from comet_ml import Experiment
from comet_ml.integration.pytorch import log_model
from accelerate import Accelerator

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.loader import DataLoader
from pepGraph_model import MixGCNBiLSTM
from pepGraph_BiLSTM import GNN_v2

import pandas as pd
import numpy as np
from tqdm import tqdm
from pdb2sql import interface

from scipy.stats import pearsonr
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
from pepGraph_utlis import seq_embedding, neighbor_search, neighbor_filter

def train_model(model, num_epochs, optimizer, train_loader, val_loader, loss_fn, accelerator, experiment, result_dir, data_log = True):
    model.reset_parameters()
    rp_train = []
    rmse_train_list = []
    best_val_loss = float('inf')

    for epoch in range(num_epochs):
        list1_train = []
        list2_train = []
        epoch_train_losses = []

        model.train()
        for graph_batch in train_loader:
            targets = graph_batch.y
            outputs = model(graph_batch)
            train_loss = loss_fn(outputs, targets)

            optimizer.zero_grad()
            accelerator.backward(train_loss)
            optimizer.step()

            epoch_train_losses.append(train_loss.item())
            targets = targets.detach().cpu()
            outputs = outputs.detach().cpu()
            list1_train=np.append(list1_train,targets)
            list2_train=np.append(list2_train,outputs)

        #accelerator.wait_for_everyone()
        if accelerator.is_main_process:
            epoch_train_losses = np.mean(epoch_train_losses)
            epoch_rp_train = np.corrcoef(list2_train, list1_train)[0,1]
            rp_train.append(np.mean(epoch_rp_train))

            x = np.array(list1_train).reshape(-1,1)
            y = np.array(list2_train).reshape(-1,1)
            epoch_train_rmse= np.sqrt(((y - x) ** 2).mean())
            rmse_train_list.append(epoch_train_rmse)

            print('Epoch  Train Loss  PCC Train  Train RMSE')
            print("{:5d}  {:10.3f}  {:9.3f}  {:10.3f}".format(
            epoch, epoch_train_losses, epoch_rp_train, epoch_train_rmse))
            
            if data_log:
                experiment.log_metric('train_loss', epoch_train_losses, step = epoch)
                experiment.log_metric('train_pcc', epoch_rp_train, step = epoch)
                experiment.log_metric('train_rmse', epoch_train_rmse, step = epoch)

        ### validation and early-stopping
        if epoch % 10 == 0:
            model.eval()
            epoch_val_losses = []
            with torch.no_grad():
                for graph_batch in val_loader:
                    targets = graph_batch.y
                    outputs = model(graph_batch)

                    all_predictions, all_targets = accelerator.gather_for_metrics((outputs, targets))
                    val_loss = loss_fn(all_predictions, all_targets)

                    epoch_val_losses.append(val_loss.item())
                val_losses_mean = np.mean(epoch_val_losses)

            accelerator.wait_for_everyone()
            if accelerator.is_main_process and data_log:
                experiment.log_metric('val_loss', val_losses_mean, step = epoch)

            if val_losses_mean < best_val_loss:
                best_val_loss = val_losses_mean
                trigger_times = 0
                accelerator.save_model(model, f'{result_dir}/best_model.pt')
            else:
                trigger_times += 1
                if trigger_times >= 3:
                    print('Early stopping')
                    break

    return rmse_train_list, rp_train

def test_model(model, test_loader):
    y_pred = []
    y_true = []
    model.eval()
    for i, graph_batch in enumerate(test_loader):
        targets = graph_batch.y.to(dtype=torch.float32)
        outputs = model(graph_batch)

        y_pred.append(outputs.cpu().detach().numpy())
        y_true.append(targets.cpu().detach().numpy())
    y_pred = np.concatenate(y_pred, axis=0)
    y_true = np.concatenate(y_true, axis=0)
    return y_true, y_pred

def main(training_args):
    ##################################### initial setting ##################################### 
    root_dir = "/home/lwang/models/HDX_LSTM/data/Fullset"
    summary_HDX_file = f'{root_dir}/merged_data_oldVer.xlsx'
    hdx_df = pd.read_excel(summary_HDX_file, sheet_name='Sheet1')
    hdx_df = hdx_df.dropna(subset=['chain_identifier'])

    pepGraph_dir = os.path.join(root_dir, 'graph_ensemble')
    result_dir = training_args['result_dir']
    if not os.path.exists(result_dir):
        os.mkdir(result_dir)

    ### data preparation ###
    apo_input = []
    complex_input = []
    hdx_df = hdx_df.drop_duplicates(subset=['apo_identifier'])
    for row_id, row in hdx_df.iterrows():
        pdb = row['apo_identifier'].strip().split('.')[0]
        pepGraph_file = f'{pepGraph_dir}/{pdb}.pt'
        if os.path.isfile(pepGraph_file):
            pepGraph_ensemble = torch.load(pepGraph_file)
            if row['complex_state'] == 'single':
                apo_input.extend(pepGraph_ensemble)
            else:
                complex_input.extend(pepGraph_ensemble)
        else:
            continue

    print('length of apo data:', len(apo_input))
    print('length of complex data:', len(complex_input))

    ### model initialization ###
    Mix_model = MixGCNBiLSTM(training_args)
    model = Mix_model

    ### training ###
    loss_fn = nn.BCELoss()    
    optimizer = torch.optim.Adam(model.parameters(), lr=config['learning_rate'], weight_decay=config['weight_decay'])

    for i in range(config['cross_validation_num']):
        print(f'-----------------------------------------')
        print(f'-         cross validation {i}          -')
        print(f'-----------------------------------------')

        train_apo, val_apo = train_test_split(apo_input, test_size=0.3, random_state=42)
        train_complex, val_complex = train_test_split(complex_input, test_size=0.3, random_state=42)

        train_set = train_apo + train_complex
        val_set = val_apo + val_complex
        train_loader = DataLoader(train_set, batch_size = config['batch_size'], shuffle=True, num_workers=config['num_workers'])
        val_loader =  DataLoader(val_set, batch_size = config['batch_size'], shuffle=False, num_workers=config['num_workers'])
        print('length of train_Set:', len(train_set))
        print('length of val_Set:', len(val_set))

        # train and save model at checkpoints
        model, optimizer, train_loader, val_loader = accelerator.prepare(model, optimizer, train_loader, val_loader)
        rmse_train_list, rp_train = train_model(
            model, config['num_epochs'], optimizer, train_loader, val_loader, loss_fn,
            accelerator, experiment, result_dir, data_log = training_args['data_log'])
        if training_args['data_log']:
            log_model(experiment, model=model, model_name = 'PEP-HDX')

    '''
        # PCC Plot
        plt.figure(figsize=(10, 6))
        plt.plot(rp_train, label='Training PCC')
        plt.title('Pearson Correlation Coefficient over epochs')
        plt.ylim(0, 1)
        plt.xlabel('Epochs')
        plt.ylabel('PCC')
        plt.legend()
        plt.savefig(f'{result_dir}/pcc_plot_{i}.png')  # Save the plot as a PNG file

        # RMSE Plot
        plt.figure(figsize=(10, 6))
        plt.plot(rmse_train_list, label='train RMSE')
        plt.title('Root Mean Square Error over epochs')
        plt.ylim(0.05, 0.2)
        plt.xlabel('Epochs')
        plt.ylabel('RMSE')
        plt.legend()
        plt.savefig(f'{result_dir}/rmse_plot_{i}.png')  # Save the plot as a PNG file

        def val_plot(y_true, y_pred, filename):
            if len(y_true) != len(y_pred):
                raise ValueError("y_true and y_pred must have the same length.")
            
            plt.figure()
            plt.scatter(y_true, y_pred, alpha=0.5, label='Data points')
            plt.xlabel('Experimental HDX')
            plt.ylabel('Predicted HDX')
            plt.xlim(0, 1)
            plt.ylim(0, 1)
            plt.plot([0, 1], [0, 1], '--', color='gray', label='Ideal correlation')
            def check_numeric(arr):
                return np.all(np.isfinite(arr))
            print(check_numeric(y_true))
            print(check_numeric(y_pred))
            #pcc, _ = pearsonr(y_true, y_pred)
            #plt.text(0.05, 0.9, f'PCC: {pcc:.2f}', transform=plt.gca().transAxes)
            #plt.savefig(f'{result_dir}/{filename}.png')

        val_apo_loader =  DataLoader(val_apo, batch_size = config['batch_size'], shuffle=False, num_workers=config['num_workers'])
        val_complex_loader =  DataLoader(val_complex, batch_size = config['batch_size'], shuffle=False, num_workers=config['num_workers'])
        val_apo_loader, val_complex_loader = accelerator.prepare(val_apo_loader, val_complex_loader)

        y_true, y_pred = test_model(model, val_apo_loader)
        val_plot(y_true, y_pred, f'val_apo_{i}')
        y_true, y_pred = test_model(model, val_complex_loader)
        val_plot(y_true, y_pred, f'val_complex_{i}')
    '''

if __name__ == "__main__":
    accelerator = Accelerator()
    device = accelerator.device
    config = {
            'num_epochs': 300,
            'batch_size': 16,
            'learning_rate': 0.001,
            'weight_decay': 5e-4,
            'GNN_type': 'GAT',
            'num_GNN_layers': 3,
            'cross_validation_num': 1,
            'num_workers': 4
    }

    training_args = {'num_hidden_channels': 10, 'num_out_channels': 20, 

            'feat_in_dim': 44, 'topo_in_dim': 42, 'num_heads': 8, 'GNN_hidden_dim': 32,
            'GNN_out_dim': 16, 'LSTM_out_dim': 16,

            'final_hidden_dim': 16,

            'drop_out': 0.5, 'num_GNN_layers': config['num_GNN_layers'], 'GNN_type': config['GNN_type'],
            'graph_hop': 'hop1', 'batch_size': config['batch_size'],
            'result_dir': '/home/lwang/models/HDX_LSTM/results/240411_BiLSTMGAT_v2_1',
            'data_log': True
    }

    os.environ["COMET_GIT_DIRECTORY"] = "/home/lwang/AI-HDX-main/ProteinComplex_HDX_prediction"  

    experiment = Experiment(
        api_key="yvWYrZuk8AhNnXgThBrgGChY4",
        project_name="CoolHdx_project",
        workspace="superchrisw"
    )
   
    if training_args['data_log']:
        experiment.log_parameters(config)

    #experiment = ''
    main(training_args)