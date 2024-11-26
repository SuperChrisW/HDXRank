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
#from pepGraph_model import MixBiLSTM_GearNet, BiLSTM, GCN

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error
from scipy.stats import spearmanr

import argparse

def train_model(model, num_epochs, optimizer, train_loader, val_loader, loss_fn, device, experiment, result_fpath, data_log, rm_feat = 'none'):
    rp_train = []
    rmse_train_list = []

    for epoch in range(num_epochs):
        list1_train = []
        list2_train = []
        epoch_train_losses = []

        model.train()
        for graph_batch in train_loader:
            graph_batch = graph_batch.to(device)
            targets = graph_batch.y
            if rm_feat == 'msa':
                node_feat = graph_batch.residue_feature[:,30:].float()
            elif rm_feat == 'sequence':
                node_feat = torch.cat([graph_batch.residue_feature[:,:30].float(), graph_batch.residue_feature[:,40:].float()], dim=1) # remove seq feat
            elif rm_feat == 'physical':
                node_feat = torch.cat([graph_batch.residue_feature[:,:40].float(), graph_batch.residue_feature[:,44:].float()], dim=1) # remove physical feat
            elif rm_feat == 'geometric':
                node_feat = torch.cat([graph_batch.residue_feature[:,:44].float(), graph_batch.residue_feature[:,56:].float()], dim=1) # remove geometric feat
            elif rm_feat == 'heteroatom':
                node_feat = graph_batch.residue_feature[:,:56].float() # remove heteroatom feat
            elif rm_feat == 'none':
                node_feat = graph_batch.residue_feature.float()
            elif rm_feat == '36feat':
                node_feat = torch.cat([graph_batch.residue_feature[:,:35].float(), graph_batch.residue_feature[:,40:41].float()], dim=1) # same as AI-HDX
            elif rm_feat == 'esm2':
                node_feat = graph_batch.residue_feature.float()
            elif rm_feat == 'hse':
                node_feat = torch.cat([graph_batch.residue_feature[:,:41].float(), graph_batch.residue_feature[:,44:].float()], dim=1)
            else:
                raise ValueError('Invalid feature type to remove')
            outputs = model(graph_batch, node_feat) # for GearNet and GearNet-Edge

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

        y = np.array(list1_train).reshape(-1,1)
        x = np.array(list2_train).reshape(-1,1)
        epoch_train_rmse= np.sqrt(((y - x) ** 2).mean())
        rmse_train_list.append(epoch_train_rmse)

        print('Epoch  Train_Loss  Train_SPR  Train RMSE')
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

                    targets = graph_batch['residue'].y
                    node_feat = graph_batch.residue_feature.float()
                    outputs = model(graph_batch, node_feat)
                
                    val_loss = loss_fn(outputs, targets)
                    epoch_val_losses.append(val_loss.item())
                val_losses_mean = np.mean(epoch_val_losses)

            if data_log:
                experiment.log_metric('val_loss', val_losses_mean, step = epoch)'''

        #if epoch % 10 == 0:
        #    save_checkpoint(model, optimizer, epoch, f'{result_fpath}_epoch{epoch}.pth')

    model.eval()
    complex_labels = []
    val_preds = []
    val_targets = []
    with torch.no_grad():
        for graph_batch in val_loader: 
            graph_batch = graph_batch.to(device)
            targets = graph_batch.y
            if rm_feat == 'msa':
                node_feat = graph_batch.residue_feature[:,30:].float()
            elif rm_feat == 'sequence':
                node_feat = torch.cat([graph_batch.residue_feature[:,:30].float(), graph_batch.residue_feature[:,40:].float()], dim=1) # remove seq feat
            elif rm_feat == 'physical':
                node_feat = torch.cat([graph_batch.residue_feature[:,:40].float(), graph_batch.residue_feature[:,44:].float()], dim=1) # remove physical feat
            elif rm_feat == 'geometric':
                node_feat = torch.cat([graph_batch.residue_feature[:,:44].float(), graph_batch.residue_feature[:,56:].float()], dim=1) # remove geometric feat
            elif rm_feat == 'heteroatom':
                node_feat = graph_batch.residue_feature[:,:56].float() # remove heteroatom feat
            elif rm_feat == 'none':
                node_feat = graph_batch.residue_feature.float()
            elif rm_feat == 'esm2':
                node_feat = graph_batch.residue_feature.float()
            elif rm_feat == 'hse':
                node_feat = torch.cat([graph_batch.residue_feature[:,:41].float(), graph_batch.residue_feature[:,44:].float()], dim=1)
            #node_feat = torch.cat([graph_batch.residue_feature[:,:35].float(), graph_batch.residue_feature[:,40:41].float()], dim=1) # same as AI-HDX
            outputs = model(graph_batch, node_feat)

            val_preds.extend(outputs.cpu().numpy())
            val_targets.extend(targets.cpu().numpy())
            complex_labels.extend(graph_batch.is_complex.cpu().numpy())
        
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

    #save_checkpoint(model, optimizer, num_epochs, f'{result_fpath}_epoch{epoch}.pth')
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

    pepGraph_dir = os.path.join(root_dir, 'graph_ensemble_simpleGearNet', training_args['cluster'])
    result_dir = training_args['result_dir']
    if not os.path.exists(result_dir):
        os.mkdir(result_dir)

    device = training_args['device']
    ### data preparation ###
    apo_input = []
    complex_input = []
    extend_apo = []
    extend_complex = []
    hdx_df = hdx_df.drop_duplicates(subset=['structure_file'])

    print('data loading')
    for row_id, row in tqdm(hdx_df.iterrows()):
        pdb = row['structure_file'].strip().split('.')[0].upper()
        pepGraph_file = f'{pepGraph_dir}/{pdb}.pt'
        extend_pepGraph = f'{pepGraph_dir}/extend_data/{pdb}.pt'
        if os.path.isfile(pepGraph_file):
            pepGraph_ensemble = torch.load(pepGraph_file)
            if row['complex_state'] == 'single':
                apo_input.extend(pepGraph_ensemble)
            else:
                complex_input.extend(pepGraph_ensemble)
        else:
            continue
        if os.path.isfile(extend_pepGraph) and training_args['use_extend']:
            extend_pepGraph_ensemble = torch.load(extend_pepGraph)
            if row['complex_state'] == 'single':
                extend_apo.extend(extend_pepGraph_ensemble)
            else:
                extend_complex.extend(extend_pepGraph_ensemble)

    if training_args['use_extend']:
        apo_input.extend(extend_apo)
        complex_input.extend(extend_complex)

    print('length of apo data:', len(apo_input))
    print('length of complex data:', len(complex_input))

    ### model initialization ###
    torch.cuda.empty_cache()
    edge_relation = {'KNNEdge':1, 'radiusEdge':1, 'seqEdge':4, 'none': 0}

    ### training ###
    for i in range(config['cross_validation_num']):
        model = GearNet(input_dim = training_args['feat_in_dim']+training_args['topo_in_dim'], hidden_dims = [512,512,512],
                            num_relation=7-edge_relation[training_args['rm_edge']], batch_norm=True, concat_hidden=True, readout='sum', activation = 'relu', short_cut=True).to(device)
        
        loss_fn = nn.BCELoss()    
        optimizer = torch.optim.Adam(model.parameters(), lr=config['learning_rate'], weight_decay=config['weight_decay'])

        print(f'-----------------------------------------')
        print(f'-         cross validation {i}            -')
        print(f'-----------------------------------------')

        train_apo, val_apo = train_test_split(apo_input, test_size=0.2, random_state=42)
        train_complex, val_complex = train_test_split(complex_input, test_size=0.2, random_state=42)


        train_set = train_apo + train_complex
        val_set = val_apo + val_complex
        train_set = data.Protein.pack(train_set+val_set)
        val_set = data.Protein.pack(val_set)

        #train_set = data.Protein.pack([item for sublist in train_apo + train_complex for item in sublist])
        #val_set = data.Protein.pack([item for sublist in val_apo + val_complex for item in sublist])
        train_set.view = 'residue'
        val_set.view = 'residue'

        train_loader = data.DataLoader(train_set, batch_size = config['batch_size'], shuffle=True, num_workers=config['num_workers'])
        val_loader =  data.DataLoader(val_set, batch_size = config['batch_size'], shuffle=False, num_workers=config['num_workers'])

        print('length of train_Set:', len(train_set))
        print('length of val_Set:', len(val_set))

        # train and save model at checkpoints
        if training_args['use_extend']:
            fname = training_args['file_name']+'_v'+str(i)+'_extend'
        else:
            fname = training_args['file_name']+'_v'+str(i)+'_'+training_args['rm_feat']
        result_fpath = os.path.join(training_args['result_dir'], fname)
        rmse_train_list, rp_train = train_model(
            model, config['num_epochs'], optimizer, train_loader, val_loader, loss_fn, device,
            experiment, result_fpath, training_args['data_log'], training_args['rm_feat'])
        if training_args['data_log']:
            log_model(experiment, model=model, model_name = 'PEP-HDX')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Train a HDX-GearNet model.')
    parser.add_argument('-cluster', type=str, required=True, default='cluster1', help='Data cluster to train the model')

    parser.add_argument('-save', type=str, default='/home/lwang/models/HDX_LSTM/results/241110_GearNet', help='Path to save the output results')
    parser.add_argument('-epoch', type=int, default= 100, help='epoch numbers to be trained')
    parser.add_argument('-cross_num', type=int, default= 1, help='cross validation numbers')
    parser.add_argument('-cuda', type=int, default= 0, help='cuda number')
    parser.add_argument('-rm_feat', type=str, default= "none", help='remove feature type from the input features')
    parser.add_argument('-rm_edge', type=str, default= "none", help='remove edge from the input features')
    args = parser.parse_args()

    cluster = args.cluster
    device = torch.device(f'cuda:{args.cuda}' if torch.cuda.is_available() else 'cpu')
    config = {
            'num_epochs':args.epoch,
            'batch_size': 16,
            'learning_rate': 0.001,
            'weight_decay': 5e-4,
            'GNN_type': f'model_GN_{cluster}',
            'num_GNN_layers': 3,
            'cross_validation_num': args.cross_num,
            'num_workers': 4,
    }

    feat_num = {"sequence": 10, "msa": 30, "physical": 4, "geometric": 12, "heteroatom": 42, 'none': 0, '36feat': 62, 'esm2': -264,
                'hse': 3}

    training_args = {'num_hidden_channels': 10, 'num_out_channels': 20, 
            'feat_in_dim': 56 - feat_num[args.rm_feat], 'topo_in_dim': 0, 'num_heads': 8, 'GNN_hidden_dim': 32,
            'GNN_out_dim': 64, 'LSTM_out_dim': 64, 'hidden_dims': [512, 512, 512],

            'final_hidden_dim': 16,

            'drop_out': 0.5, 'num_GNN_layers': config['num_GNN_layers'], 'GNN_type': config['GNN_type'],
            'graph_hop': 'hop1', 'batch_size': config['batch_size'],
            'result_dir': args.save,
            'file_name': f'model_GN56_{cluster}',
            'data_log': False,
            'device': device,
            'cluster': cluster,
            'rm_feat': args.rm_feat,
            'rm_edge': args.rm_edge,
            'use_extend': False, # use extra 2671 synthetic peptides with RFU in range [0.6, 1.0]
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