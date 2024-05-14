### training ###
import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

from comet_ml import Experiment
from comet_ml.integration.pytorch import log_model

import torch
import torch.nn as nn
import torch.nn.functional as F
from torchdrug import data
from GearNet import GearNet
from pepGraph_model import MixBiLSTM_GearNet

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split


def train_model(model, num_epochs, optimizer, train_loader, val_loader, loss_fn, device, experiment, result_fpath, data_log = True):
    #model.reset_parameters()
    rp_train = []
    rmse_train_list = []
    best_val_loss = float('inf')

    for epoch in range(num_epochs):
        list1_train = []
        list2_train = []
        epoch_train_losses = []

        model.train()
        for graph_batch in train_loader:
            graph_batch = graph_batch.to(device)
            targets = graph_batch.y
            outputs = model(graph_batch, graph_batch.residue_feature.float())
            #outputs = model(graph_batch)

            train_loss = loss_fn(outputs, targets)
            optimizer.zero_grad()
            train_loss.backward()
            optimizer.step()

            epoch_train_losses.append(train_loss.item())
            targets = targets.detach().cpu()
            outputs = outputs.detach().cpu()
            list1_train=np.append(list1_train,targets)
            list2_train=np.append(list2_train,outputs)

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
                    graph_batch = graph_batch.to(device)
                    targets = graph_batch.y
                    outputs = model(graph_batch, graph_batch.residue_feature.float())
                    #outputs = model(graph_batch)

                    val_loss = loss_fn(outputs, targets)
                    epoch_val_losses.append(val_loss.item())
                val_losses_mean = np.mean(epoch_val_losses)

            if data_log:
                experiment.log_metric('val_loss', val_losses_mean, step = epoch)

            if val_losses_mean < best_val_loss:
                best_val_loss = val_losses_mean
                trigger_times = 0
                torch.save(model.state_dict(), f'{result_fpath}.pth')
            else:
                trigger_times += 1
                if trigger_times >= 5:
                    print('Early stopping')
                    break

    return rmse_train_list, rp_train

def main(training_args):
    ##################################### initial setting ##################################### 
    root_dir = "/home/lwang/models/HDX_LSTM/data/Fullset"
    summary_HDX_file = f'{root_dir}/merged_data_oldVer.xlsx'
    hdx_df = pd.read_excel(summary_HDX_file, sheet_name='Sheet1')
    hdx_df = hdx_df.dropna(subset=['chain_identifier'])

    pepGraph_dir = os.path.join(root_dir, 'graph_ensemble_GearNetEdge_minus')
    result_dir = training_args['result_dir']
    result_fpath = os.path.join(training_args['result_dir'], training_args['file_name'])
    if not os.path.exists(result_dir):
        os.mkdir(result_dir)

    device = training_args['device']
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
    torch.cuda.empty_cache()

    #GearNet
    model = GearNet(input_dim = training_args['feat_in_dim']+training_args['topo_in_dim'], hidden_dims = [512,512,512],
                    num_relation=7, batch_norm=True, concat_hidden=True, readout='sum', activation = 'relu', short_cut=True).to(device)
    
    #GearNet-Edge
    #model = GearNet(input_dim=training_args['feat_in_dim']+training_args['topo_in_dim'], hidden_dims=[512, 512, 512], 
    #                          num_relation=7, edge_input_dim=59, num_angle_bin=8,
    #                          batch_norm=True, concat_hidden=True, short_cut=True, readout="sum", activation = 'relu').to(device)

    #MixBiLSTM_GearNet
    #model = MixBiLSTM_GearNet(training_args).to(device)

    ### training ###
    loss_fn = nn.BCELoss()    
    optimizer = torch.optim.Adam(model.parameters(), lr=config['learning_rate'], weight_decay=config['weight_decay'])

    for i in range(config['cross_validation_num']):
        print(f'-----------------------------------------')
        print(f'-         cross validation {i}          -')
        print(f'-----------------------------------------')

        train_apo, val_apo = train_test_split(apo_input, test_size=0.3, random_state=42)
        train_complex, val_complex = train_test_split(complex_input, test_size=0.3, random_state=42)

        train_set = data.Protein.pack(train_apo + train_complex)
        val_set = data.Protein.pack(val_apo + val_complex)
        train_set.view = 'residue'
        val_set.view = 'residue'

        train_loader = data.DataLoader(train_set, batch_size = config['batch_size'], shuffle=True, num_workers=config['num_workers'])
        val_loader =  data.DataLoader(val_set, batch_size = config['batch_size'], shuffle=False, num_workers=config['num_workers'])
        print('length of train_Set:', len(train_set))
        print('length of val_Set:', len(val_set))

        # train and save model at checkpoints
        rmse_train_list, rp_train = train_model(
            model, config['num_epochs'], optimizer, train_loader, val_loader, loss_fn, device,
            experiment, result_fpath, data_log = training_args['data_log'])
        if training_args['data_log']:
            log_model(experiment, model=model, model_name = 'PEP-HDX')

if __name__ == "__main__":

    device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')
    config = {
            'num_epochs':300,
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
            'result_dir': '/home/lwang/models/HDX_LSTM/results/240509_avgRFU',
            'file_name': 'best_model_GearNet_minus',
            'data_log': True,
            'device': device
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
    torch.cuda.empty_cache()