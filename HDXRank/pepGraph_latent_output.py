### prediction ###
import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

import torch
import torch.nn as nn
import torch.nn.functional as F
from torchdrug import data
from torch_geometric.loader import DataLoader as tg_DataLoader

import pandas as pd
import numpy as np
from GearNet import GearNet
from tqdm import tqdm

def test_model(model, test_loader, device, graph_type = 'GearNetEdge'):
    y_pred = []
    y_true = []
    range_list = []
    chain_list = []
    graph_embedding_matrix = []
    model.eval()
    try:
        with torch.no_grad():
            for i, graph_batch in enumerate(test_loader):
                graph_batch = graph_batch.to(device)
                targets = graph_batch.y
                node_feat = graph_batch.residue_feature.float()
                outputs, graph_embed = model(graph_batch, node_feat)

                graph_embedding_matrix.append(graph_embed.cpu().detach().numpy())
                range_list.extend(graph_batch.range.cpu().detach().numpy())
                chain_list.extend(graph_batch.chain.cpu().detach().numpy())
                y_pred.append(outputs.cpu().detach().numpy())
                y_true.append(targets.cpu().detach().numpy())
            graph_embedding_matrix = np.concatenate(graph_embedding_matrix, axis=0)
            print('graph_embedding_matrix:', graph_embedding_matrix.shape)
            y_pred = np.concatenate(y_pred, axis=0) if len(y_pred) > 0 else []
            y_true = np.concatenate(y_true, axis=0) if len(y_true) > 0 else []
        return y_true, y_pred, range_list, chain_list, graph_embedding_matrix
    except Exception as e:
        print(e)
        return None, None, None, None

def run_prediction(model, hdx_df, merged_config, pepGraph_dir, device):
    print('data loading')
    input_data = []
    for row_id, row in tqdm(hdx_df.iterrows()):
        pdb = row['structure_file'].strip().split('.')[0].upper()
        pepGraph_file = f'{pepGraph_dir}/{pdb}.pt'
        if os.path.isfile(pepGraph_file):
            print('loading:', pdb)
            pepGraph_ensemble = torch.load(pepGraph_file)
            input_data.extend(pepGraph_ensemble)
        else:
            continue
    pred_set = data.Protein.pack(input_data)
    pred_set.view = 'residue'
    test_loader = data.DataLoader(pred_set, batch_size = merged_config['batch_size'], shuffle=False, num_workers=merged_config['num_workers'])
    y_true, y_pred, range_list, chain_list, graph_embedding_matrix = test_model(model, test_loader, device, merged_config['graph_type'])
    return y_true, y_pred, range_list, chain_list, graph_embedding_matrix

def load_model(model_path, config, accelerator):
    model = GearNet(config)
    unwrapped_model = accelerator.unwrap_model(model)
    path_to_checkpoint = os.path.join(model_path, "pytorch_model.bin")
    unwrapped_model.load_state_dict(torch.load(path_to_checkpoint))
    return unwrapped_model

def main(save_args):
    ##################################### initial setting #####################################
    root_dir = "/home/lwang/models/HDX_LSTM/data/Latest_set"
    summary_HDX_file = f'{root_dir}/merged_data-NAcomplex.xlsx'
    hdx_df = pd.read_excel(summary_HDX_file, sheet_name='Sheet1')
    hdx_df = hdx_df.dropna(subset=['structure_file'])

    pepGraph_dir = os.path.join(root_dir, 'graph_ensemble_simpleGearNet', save_args['cluster'])
    model_path = save_args['model_path']
    device = save_args['device']
    ##################################### config setting #####################################
    config = {
    #training setting
    'num_epochs': 8, # epochs for finetuning
    'batch_size': 16,
    'learning_rate': 0.001,
    'weight_decay': 5e-4,
    'GNN_type': 'GearNet',
    'num_GNN_layers': 3,
    'cross_validation_num': 1,
    'num_workers': 4,
    'data_log': False,  

    #model setting
    'num_hidden_channels': 10,
    'num_out_channels': 20,
    'feat_in_dim': 56, #bilstm 36: 56-62, GVP 1280: 1280, GearNet: 56
    'topo_in_dim': 0,
    'num_heads': 8,
    'GNN_hidden_dim': 32,
    'GNN_out_dim': 64,
    'LSTM_out_dim': 64,
    'final_hidden_dim': 16,
    'hidden_dims': [512, 512, 512],
    'drop_out': 0.5,
    'graph_hop': 'hop1'
    }
    merged_config = {**config, **save_args}

    ##################################### model setting #####################################
    #GearNet
    model = GearNet(input_dim = merged_config['feat_in_dim']+merged_config['topo_in_dim'], hidden_dims = [512,512,512],
                    num_relation=7, batch_norm=True, concat_hidden=True, readout='sum', activation = 'relu', short_cut=True)
 
    model_state_dict = torch.load(model_path, map_location=device)
    model_state_dict = model_state_dict['model_state_dict'] # add if saved as checkpoint
    model.load_state_dict(model_state_dict)
    model = model.to(device)
    print('model loaded successfully!')

    # run_prediciton here 
    merged_config['load_proteins'] = merged_config['prediction_protein']
    y_true, y_pred, range_, chain, graph_embedding_matrix = run_prediction(model, hdx_df, merged_config, pepGraph_dir, device)
    print('final graph_embed_list shape:', graph_embedding_matrix.shape)

    data = {
    'Y_True': np.array(y_true),
    'Y_Pred': np.array(y_pred),
    'Graph_Embedding': graph_embedding_matrix
    }
    np.save(f'{merged_config["result_dir"]}/graph_embedding.npy', data)

if __name__ == "__main__":
    root_dir = "/home/lwang/models/HDX_LSTM/data/Latest_set"
    summary_HDX_file = f'{root_dir}/merged_data-NAcomplex.xlsx'
    hdx_df = pd.read_excel(summary_HDX_file, sheet_name='Sheet1')
    hdx_df = hdx_df.dropna(subset=['structure_file'])

    cluster = 'cluster1_8A_manual_rescale'
    graph_type = 'GearNet'
    epoch = [100]
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

    protein_list = list(hdx_df['structure_file'].unique())
    save_args = {
        #save setting
        'plot_title': 'test_set',
        'date': '240730',
        'model_name': 'GearNet',
        'version': 'v1',
        'result_dir': f'/home/lwang/models/HDX_LSTM/data/hdock/prediction/',
        'prediction_csv': None,
        'cluster':cluster,
        'graph_type': graph_type,  ## 'GearNetEdge' or 'GVP' or 'Bilstm36'
        'device': device,

        #prediction setting
        'prediction_protein': protein_list,
        'model_path': None,
        #plot setting
        'slide_window': 1,
        'show_truth': True,
        'show_pred': True,
        }
    
    for i in range(len(epoch)):
        #240920 GearNet better performance in test set
        save_args['model_path'] = f'/home/lwang/models/HDX_LSTM/results/240918_GVP/model_GN56_cluster1_8A_manual_rescale_v0_epoch{epoch[i]-1}.pth'
        save_args['result_dir'] = f'/home/lwang/models/HDX_LSTM/results/240918_GVP'
        if not os.path.exists(save_args['result_dir']):
            os.mkdir(save_args['result_dir'])
        main(save_args)
    