### training ###
import os
import sys
from tqdm import tqdm

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
from pepGraph_model import MixBiLSTM_GearNet, BiLSTM

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split


def train_model(model, num_epochs, optimizer, train_loader, val_loader, loss_fn, device, experiment, result_fpath, data_log):
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
            new_embedding = torch.cat((graph_batch.seq_embedding[:,:,:35],graph_batch.seq_embedding[:,:,39:40]), dim=-1)
            outputs = model(new_embedding.unsqueeze(1)) # for MixBiLSTM_GearNet
            #outputs = model(graph_batch.seq_embedding.unsqueeze(1)) # for BiLSTM

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

    torch.save(model.state_dict(), f'{result_fpath}.pth')
    return rmse_train_list, rp_train

def main(training_args):
    ##################################### initial setting ##################################### 
    root_dir = "/home/lwang/models/HDX_LSTM/data/Latest_set"
    summary_HDX_file = f'{root_dir}/merged_data.xlsx'
    hdx_df = pd.read_excel(summary_HDX_file, sheet_name='Sheet1')
    hdx_df = hdx_df.dropna(subset=['structure_file'])

    pepGraph_dir = os.path.join(root_dir, 'graph_ensemble_GearNetEdge', training_args['cluster'])
    result_dir = training_args['result_dir']
    result_fpath = os.path.join(training_args['result_dir'], training_args['file_name'])
    if not os.path.exists(result_dir):
        os.mkdir(result_dir)

    device = training_args['device']
    ### data preparation ###
    hdx_df = hdx_df.drop_duplicates(subset=['structure_file'])

    apo_input = []
    complex_input = []
    print('data loading')
    for row_id, row in tqdm(hdx_df.iterrows()):
        pdb = row['structure_file'].strip().split('.')[0].upper()
        pepGraph_file = f'{pepGraph_dir}/{pdb}.pt'
        if os.path.isfile(pepGraph_file):
            pepGraph_ensemble = torch.load(pepGraph_file)
            if pepGraph_ensemble[0].seq_embedding.shape[1] != 98:
                continue
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

    #MixBiLSTM_GearNet
    #model = MixBiLSTM_GearNet(training_args).to(device)

    #BiLSTM
    training_args['feat_in_dim'] = 36
    training_args['topo_in_dim'] = 0
    model = BiLSTM(training_args).to(device)

    ### training ###
    for i in range(config['cross_validation_num']):
        model.reset_parameters()
        loss_fn = nn.BCELoss()    
        optimizer = torch.optim.Adam(model.parameters(), lr=config['learning_rate'], weight_decay=config['weight_decay'])

        print(f'-----------------------------------------')
        print(f'-         cross validation {i}          -')
        print(f'-----------------------------------------')

        train_apo, val_apo = train_test_split(apo_input, test_size=0.2, random_state=42)
        train_complex, val_complex = train_test_split(complex_input, test_size=0.2, random_state=42)

        train_set = data.Protein.pack(train_apo + train_complex)
        val_set = data.Protein.pack(val_apo + val_complex)
        train_set.view = 'residue'
        val_set.view = 'residue'

        #train_set = torch.load(f'{pepGraph_dir}/train_lstm.pt')
        #val_set = torch.load(f'{pepGraph_dir}/val.pt')

        train_loader = data.DataLoader(train_set, batch_size = config['batch_size'], shuffle=True, num_workers=config['num_workers'])
        val_loader =  data.DataLoader(val_set, batch_size = config['batch_size'], shuffle=False, num_workers=config['num_workers'])
    
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
    cluster = 'cluster2'
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    config = {
            'num_epochs':60,
            'batch_size': 16,
            'learning_rate': 0.001,
            'weight_decay': 5e-4,
            'GNN_type': f'model_Bilstm_epoch60_{cluster}',
            'num_GNN_layers': 3,
            'cross_validation_num': 5,
            'num_workers': 4,
    }

    feat_num = {"sequence": 10, "msa": 30, "physical": 4, "geometric": 12, "heteroatom": 42, 'none': 0}

    training_args = {'num_hidden_channels': 10, 'num_out_channels': 20, 
            'feat_in_dim': 56 - feat_num["none"], 'topo_in_dim': 42, 'num_heads': 8, 'GNN_hidden_dim': 32,
            'GNN_out_dim': 64, 'LSTM_out_dim': 64,

            'final_hidden_dim': 16,

            'drop_out': 0.5, 'num_GNN_layers': config['num_GNN_layers'], 'GNN_type': config['GNN_type'],
            'graph_hop': 'hop1', 'batch_size': config['batch_size'],
            'result_dir': '/home/lwang/models/HDX_LSTM/results/240619_GearNetEdge',
            'file_name': f'model_Bilstm36_epoch60_{cluster}',
            'data_log': True,
            'device': device,
            'cluster': cluster
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