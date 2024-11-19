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
from torch_geometric.loader import DataLoader
from pepGraph_BiLSTM import GAT

import pandas as pd
import numpy as np
from tqdm import tqdm
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error
from scipy.stats import spearmanr


def train_model(model, num_epochs, optimizer, train_loader, val_loader, loss_fn, device, experiment, result_fpath, data_log = True):
    rp_train = []
    rmse_train_list = []

    for epoch in range(num_epochs):
        list1_train = []
        list2_train = []
        epoch_train_losses = []

        model.train()
        for graph_batch in train_loader:
            optimizer.zero_grad()
            graph_batch = graph_batch.to(device)
            targets = graph_batch['residue'].y.to(torch.float32)
            graph_batch['residue'].node_s['residue'] = torch.cat((graph_batch['residue'].node_s['residue'][:,:35],graph_batch['residue'].node_s['residue'][:,40:41]), dim=-1)
            outputs = model(graph_batch)

            train_loss = loss_fn(outputs, targets)
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
        '''if epoch % 5 == 0:
            model.eval()
            epoch_val_losses = []
            with torch.no_grad():
                for graph_batch in val_loader:
                    graph_batch = graph_batch.to(device)
                    targets = graph_batch['residue'].y.to(torch.float32)
                    outputs = model(graph_batch)
                        
                    val_loss = loss_fn(outputs, targets)
                    epoch_val_losses.append(val_loss.item())
                val_losses_mean = np.mean(epoch_val_losses)

            if data_log:
                experiment.log_metric('val_loss', val_losses_mean, step = epoch)

        if epoch % 10 == 0:
            save_checkpoint(model, optimizer, epoch, f'{result_fpath}_epoch{epoch}.pth')'''

    #torch.save(model.state_dict(), f'{result_fpath}.pth')
    model.eval()
    complex_labels = []
    val_preds = []
    val_targets = []
    with torch.no_grad():
        for batch in val_loader: 
            batch = batch.to(device)
            targets = batch['residue'].y.to(torch.float32)
            batch['residue'].node_s['residue'] = torch.cat((batch['residue'].node_s['residue'][:,:35],batch['residue'].node_s['residue'][:,40:41]), dim=-1)
            outputs = model(batch)
            val_preds.extend(outputs.cpu().numpy())
            val_targets.extend(targets.cpu().numpy())
            complex_labels.extend(batch['residue'].is_complex.cpu().numpy())
        
    val_preds = np.array(val_preds)
    val_targets = np.array(val_targets)
    complex_labels = np.array(complex_labels)
    single_mask = (complex_labels == 0)
    complex_mask = (complex_labels == 1)
    single_rmse = np.sqrt(mean_squared_error(val_preds[single_mask], val_targets[single_mask]))
    complex_rmse = np.sqrt(mean_squared_error(val_preds[complex_mask], val_targets[complex_mask]))
    total_rmse = np.sqrt(mean_squared_error(val_preds, val_targets))
    single_spr = spearmanr(val_preds[single_mask], val_targets[single_mask])[0]
    complex_spr = spearmanr(val_preds[complex_mask], val_targets[complex_mask])[0]
    total_spr = spearmanr(val_preds, val_targets)[0]
    single_r2 = r2_score(val_preds[single_mask], val_targets[single_mask])
    complex_r2 = r2_score(val_preds[complex_mask], val_targets[complex_mask])
    total_r2 = r2_score(val_preds, val_targets)
    print('------------ Model Validation ------------')
    print('Single RMSE   Single SPR   Single R2')
    print("{:10.3f}  {:9.3f}  {:10.3f}".format(single_rmse, single_spr, single_r2))
    print('Complex RMSE  Complex SPR  Complex R2')
    print("{:10.3f}  {:9.3f}  {:10.3f}".format(complex_rmse, complex_spr, complex_r2))
    print('Total RMSE    Total SPR    Total R2')
    print("{:10.3f}  {:9.3f}  {:10.3f}".format(total_rmse, total_spr, total_r2))

    save_checkpoint(model, optimizer, num_epochs, f'{result_fpath}_epoch{epoch}.pth')
    return rmse_train_list, rp_train

def save_checkpoint(model, optimizer, epoch, file_path):
    checkpoint = {
        'epoch': epoch,
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict()
    }
    torch.save(checkpoint, file_path)

def main(training_args):
   ##################################### initial setting ##################################### 
    root_dir = "/home/lwang/models/HDX_LSTM/data/Latest_set"
    summary_HDX_file = f'{root_dir}/merged_data-NAcomplex.xlsx'
    hdx_df = pd.read_excel(summary_HDX_file, sheet_name='Sheet1')
    hdx_df = hdx_df.dropna(subset=['structure_file'])

    pepGraph_dir = os.path.join(root_dir, 'graph_ensemble_simpleGVP', training_args['cluster'])
    result_dir = training_args['result_dir']
    if not os.path.exists(result_dir):
        os.mkdir(result_dir)

    device = training_args['device']
    ### data preparation ###
    apo_input = []
    complex_input = []
    hdx_df = hdx_df.drop_duplicates(subset=['structure_file'])

    print('data loading')
    for row_id, row in tqdm(hdx_df.iterrows()):
        pdb = row['structure_file'].strip().split('.')[0].upper()
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

    ### training ###
    loss_fn = nn.BCELoss()
    for i in range(config['cross_validation_num']):
        model = GAT(training_args).to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=config['learning_rate'], weight_decay=config['weight_decay'])

        print(f'-----------------------------------------')
        print(f'-         cross validation {i}            -')
        print(f'-----------------------------------------')

        train_apo, val_apo = train_test_split(apo_input, test_size=0.2, random_state=42)
        train_complex, val_complex = train_test_split(complex_input, test_size=0.2, random_state=42)
        print('length of train_apo:', len(train_apo))
        print('length of val_apo:', len(val_apo))
        print('length of train_complex:', len(train_complex))
        print('length of val_complex:', len(val_complex))

        train_set = train_apo + train_complex
        val_set = val_apo + val_complex

        train_loader = DataLoader(train_set, batch_size = config['batch_size'], shuffle=True, num_workers=config['num_workers'])
        val_loader =  DataLoader(val_set, batch_size = config['batch_size'], shuffle=False, num_workers=config['num_workers'])
        print('length of train_Set:', len(train_set))
        print('length of val_Set:', len(val_set))

        # train and save model at checkpoints
        fname = training_args['file_name']+'_v'+str(i)
        result_fpath = os.path.join(training_args['result_dir'], fname)
        rmse_train_list, rp_train = train_model(
            model, config['num_epochs'], optimizer, train_loader, val_loader, loss_fn, device,
            experiment, result_fpath, training_args['data_log'])
        if training_args['data_log']:
            log_model(experiment, model=model, model_name = 'PEP-HDX')

if __name__ == "__main__":
    cluster = 'cluster1_8A_manual_rescale'
    model_name = 'GAT36'
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    config = {
            'num_epochs':100,
            'batch_size': 16,
            'learning_rate': 0.001,
            'weight_decay': 5e-4,
            'GNN_type': f'model_{model_name}_{cluster}_radiusEdge+self',
            'num_GNN_layers': 3,
            'cross_validation_num': 1,
            'num_workers': 4,
    }

    training_args = {'num_hidden_channels': 10, 'num_out_channels': 20, 
            'feat_in_dim': 36, 'topo_in_dim': 0, 'num_heads': 8, 'GNN_hidden_dim': 32,
            'GNN_out_dim': 16, 'LSTM_out_dim': 16,

            'final_hidden_dim': 16,

            'drop_out': 0.5, 'num_GNN_layers': config['num_GNN_layers'], 'GNN_type': config['GNN_type'],
            'graph_hop': 'hop1', 'batch_size': config['batch_size'],
            'result_dir': '/home/lwang/models/HDX_LSTM/results/240918_GVP',
            'file_name': f'model_{model_name}_{cluster}_radiusEdge+self',
            'data_log': False,
            'device': device,
            'cluster': cluster
    }

    '''os.environ["COMET_GIT_DIRECTORY"] = "/home/lwang/AI-HDX-main/ProteinComplex_HDX_prediction"  

    experiment = Experiment(
        api_key="yvWYrZuk8AhNnXgThBrgGChY4",
        project_name="LuckHdx_project",
        workspace="superchrisw"
    )
   
    if training_args['data_log']:
        experiment.log_parameters(config)'''

    experiment = ''
    main(training_args)
    torch.cuda.empty_cache()