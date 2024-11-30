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
import itertools
from GearNet import GearNet
from pepGraph_model import GCN, BiLSTM
from GVP_model import MQAModel

from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt
from sklearn.metrics import mean_absolute_error, mean_squared_error
from math import sqrt
from tqdm import tqdm

def test_model(model, test_loader, device, graph_type = 'GearNetEdge'):
    y_pred = []
    y_true = []
    range_list = []
    chain_list = []
    model.eval()
    try:
        for i, graph_batch in enumerate(test_loader):
                graph_batch = graph_batch.to(device)
                if graph_type == 'GearNet':
                    targets = graph_batch.y
                    #node_feat = graph_batch.residue_feature[:,:56].float() # remove heteroatom encoding
                    #remove hse encoding
                    node_feat = torch.cat([graph_batch.residue_feature[:,:41].float(), graph_batch.residue_feature[:,44:].float()], dim=1)

                    #node_feat = graph_batch.residue_feature.float()
                    outputs = model(graph_batch, node_feat)
                    range_list.extend(graph_batch.range.cpu().detach().numpy())
                    chain_list.extend(graph_batch.chain.cpu().detach().numpy())
                    y_pred.append(outputs.cpu().detach().numpy())
                    y_true.append(targets.cpu().detach().numpy())
                elif graph_type == 'GVP':
                    targets = graph_batch['residue'].y
                    outputs = model((graph_batch['residue'].node_s, graph_batch['residue'].node_v), graph_batch.edge_index_dict, 
                                    (graph_batch['residue'].edge_s, graph_batch['residue'].edge_v), batch = graph_batch['residue'].batch)
                    range_list.extend(graph_batch['residue'].range.cpu().tolist())
                    chain_list.extend(graph_batch['residue'].chain.cpu().tolist())
                    y_pred.append(outputs.cpu().detach().numpy())
                    y_true.append(targets.cpu().detach().numpy())
                elif graph_type == 'Bilstm36': # load GVP data
                    targets = graph_batch['residue'].y
                    new_embedding = graph_batch['residue'].seq_embedding[:,:36]
                    outputs = model(new_embedding)
                    range_list.extend(graph_batch['residue'].range.cpu().detach().numpy())
                    chain_list.extend(graph_batch['residue'].chain.cpu().detach().numpy())
                    y_pred.append(outputs.cpu().detach().numpy())
                    y_true.append(targets.cpu().detach().numpy())
                else:
                    raise ValueError('Invalid graph type')

        y_pred = np.concatenate(y_pred, axis=0) if len(y_pred) > 0 else []
        y_true = np.concatenate(y_true, axis=0) if len(y_true) > 0 else []
        return y_true, y_pred, range_list, chain_list
    except Exception as e:
        print(e)
        return None, None, None, None

def plot_result(y_true, y_pred, range_list, batch_list, chain_list, args):
    date = args['date']
    model_name = args['model_name']
    version = args['version']
    dataset = args['plot_title']
    save_dir = args['result_dir']
    slide_window = args['slide_window']

    ### evaluation ###
    plt.figure(figsize=(20, 10))
    colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'])
    markers = itertools.cycle(['v', 'o', ',', '.', '^', '<', '>', '1', '2', '3', '4', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd', '|', '_'])
    def slide_mean(data, window_size=5):
        return np.convolve(data, np.ones(window_size)/window_size, mode='same') 

    target_batch = np.unique(batch_list)
    for i, batch in enumerate(target_batch):
        mask = (np.array(batch_list) == batch)
        
        range_list = np.array(range_list).reshape(-1, 2)
        pep_center = np.mean(range_list, axis=1)
        chain_list = np.array(chain_list)
        x_strings = np.array([f'{chain}_{int(start)}-{int(end)}' for i, (start, end) in enumerate(range_list) for chain in chain_list[i]])

        temp_y_true = y_true.flatten()[mask]
        temp_y_pred = y_pred.flatten()[mask]
        temp_pep_center = pep_center[mask]
        temp_x_strings = x_strings[mask]

        print(f'{batch} prediction:')
        pcc = pearsonr(temp_y_true, temp_y_pred)
        mae = mean_absolute_error(temp_y_true, temp_y_pred)
        rmse = sqrt(mean_squared_error(temp_y_true, temp_y_pred))
        spR = spearmanr(temp_y_true, temp_y_pred)
        print('Mean Absolute Error (MAE):', mae)
        print('PCC:', pcc[0])
        print('Root Mean Squared Error (RMSE):', rmse)
        print('Spearman:', spR[0])

        temp_pep_center, temp_y_true, temp_y_pred, temp_x_strings = zip(*sorted(zip(temp_pep_center, temp_y_true, temp_y_pred, temp_x_strings)))
        temp_y_true = slide_mean(temp_y_true, window_size=slide_window)
        temp_y_pred = slide_mean(temp_y_pred, window_size=slide_window)
        #temp_x_label = slide_mean(temp_x_label, window_size=slide_window)

        temp_x_strings = temp_x_strings[1:-1]
        temp_pep_center = temp_pep_center[1:-1]
        temp_y_true = temp_y_true[1:-1]
        temp_y_pred = temp_y_pred[1:-1]
        x_index = [i * 1 for i in range(len(temp_pep_center))]
        if args['show_pred']:
            plt.plot(x_index, temp_y_pred, label=f'{target_batch[i]}', color=next(colors), marker=next(markers), 
                    linestyle='--', linewidth=3, markersize=5, alpha = 1)
        if args['show_truth']:
            plt.plot(x_index, temp_y_true, label=f'{target_batch[i]}_truth', color='m', marker='*', linestyle='-', 
                    linewidth=3, markersize=5, alpha = 1)

    plt.xticks(x_index, temp_x_strings, rotation=90, ha='right')
    #plt.axhline(0, color='black', linewidth=0.5)
    plt.title(f'Comparison of True and Predicted HDX: {dataset}')  # Adding chart title
    plt.ylabel('HDX level')  # Labeling the y-axis
    plt.legend(loc='best')
    plt.grid(True) 
    plt.savefig(f'{save_dir}/{date}_{model_name}_{version}_{dataset}.png')  

def run_prediction(model, merged_config, pepGraph_dir, device):
    ##################################### data preparation #####################################
    test_loader = load_data(pepGraph_dir, merged_config)
    #try:
    y_true, y_pred, range_list, chain_list = test_model(model, test_loader, device, merged_config['graph_type'])
    return y_true, y_pred, range_list, chain_list
    #except Exception as e:
    #    print(e)
    #    return None, None, None, None

def load_model(model_path, config, accelerator):
    model = GearNet(config)
    unwrapped_model = accelerator.unwrap_model(model)
    path_to_checkpoint = os.path.join(model_path, "pytorch_model.bin")
    unwrapped_model.load_state_dict(torch.load(path_to_checkpoint))
    return unwrapped_model

def load_data(pepGraph_dir, merged_config):
    input_data = []
    load_proteins = merged_config['load_proteins'] if isinstance(merged_config['load_proteins'], list) else [merged_config['load_proteins']]
    for protein in load_proteins:
        pdb = protein.upper()
        pepGraph_file = os.path.join(pepGraph_dir, f'{pdb}.pt')
        if os.path.isfile(pepGraph_file):
            pepGraph_ensemble = torch.load(pepGraph_file) # list of graphs
            input_data.extend(pepGraph_ensemble)
        else:
            print('file not exist:', pepGraph_file)
    if merged_config['graph_type'] == 'GearNet':
        test_set = data.Protein.pack(input_data)
        test_set.view = 'residue'
        test_loader = data.DataLoader(test_set, batch_size = merged_config['batch_size'], shuffle=False, num_workers=merged_config['num_workers'])
    elif merged_config['graph_type'] == 'GVP' or merged_config['graph_type'] == 'Bilstm36':
        test_set = input_data
        test_loader = tg_DataLoader(input_data, batch_size=merged_config['batch_size'], shuffle=False, num_workers=merged_config['num_workers'])
    else:
        raise ValueError('Invalid graph type')
    return test_loader

def main(save_args):
    ##################################### initial setting #####################################
    root_dir = "/home/lwang/models/HDX_LSTM/data/COVID_SPIKE"
    #pepGraph_dir = os.path.join(root_dir, f"graph_ensemble_{save_args['graph_type']}", f"{save_args['protein_name']}", f"{save_args['cluster']}")
    pepGraph_dir = os.path.join(root_dir, f"graph_ensemble_GearNet", f"{save_args['protein_name']}", f"{save_args['cluster']}")
    model_path = save_args['model_path']

    device = save_args['device']
    ##################################### config setting #####################################
    config = {
    #training setting
    'num_epochs': 8, # epochs for finetuning
    'batch_size': 128,
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
    'feat_in_dim': 53, #bilstm 36: 56-62, GVP 1280: 1280, GearNet: 56 , GearNet-hse: 53
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
    
    #GearNet-Edge
    #model = GearNet(input_dim=merged_config['feat_in_dim']+merged_config['topo_in_dim'], hidden_dims=[512, 512, 512], 
    #                          num_relation=7, edge_input_dim=59, num_angle_bin=8,
    #                          batch_norm=True, concat_hidden=True, short_cut=True, readout="sum", activation = 'relu').to(device)

    #BiLSTM
    #model = BiLSTM(merged_config).to(device)

    # GCN
    #model = GCN(merged_config).to(device)

    # GVP-GNN
    '''node_in_dim = (56,3) # 78 or 98 or 1280
    node_h_dim = (1280, 12) # 392 or 1280
    edge_in_dim = (32,1)
    edge_h_dim = (128, 4)
    #metadata = (['residue'], [('residue', 'radius_edge', 'residue')])
    metadata = (['residue'], [('residue', 'radius_edge', 'residue'), ('residue', 'knn_edge', 'residue'), ('residue', 'forward_1_edge', 'residue'), ('residue', 'backward_1_edge', 'residue'),
                               ('residue', 'forward_2_edge', 'residue'), ('residue', 'backward_2_edge', 'residue'), ('residue', 'self_edge', 'residue')])
    model = MQAModel(node_in_dim, node_h_dim, 
        edge_in_dim, edge_h_dim, metadata,
        num_layers=config['num_GNN_layers'], drop_rate=config['drop_out'])'''
    
    model_state_dict = torch.load(model_path, map_location=device)
    model_state_dict = model_state_dict['model_state_dict'] # add if saved as checkpoint
    model.load_state_dict(model_state_dict)
    model = model.to(device)
    print('model loaded successfully!')

    # run_prediciton here 
    batch_list = []
    chain_list = []
    range_list = []
    y_true_list = []
    y_pred_list = []

    for protein in merged_config['prediction_protein']:
        print(f'Predicting {protein}...')
        merged_config['load_proteins'] = protein
        y_true, y_pred, range_, chain = run_prediction(model, merged_config, pepGraph_dir, device)
        if y_true is None:
            print(f'Error in predicting {protein}...')
            continue

        batch_list.extend([protein] * len(y_true))
        chain_list.extend(chain)
        range_list.extend(range_)
        y_true_list.extend(y_true)
        y_pred_list.extend(y_pred)

        '''
        pearsonr_ = pearsonr(y_true, y_pred)
        spearmanr_ = spearmanr(y_true, y_pred)
        rmse = sqrt(mean_squared_error(y_true, y_pred))
        mae = mean_absolute_error(y_true, y_pred)
        print('Pearson:', pearsonr_)
        print('Spearman:', spearmanr_)
        print('RMSE:', rmse)
        print('MAE:', mae)'''

    range_list = np.array(range_list).reshape(-1, 2)
    chain_list = np.array(chain_list)
    x_strings = np.array([f'{int(start)}-{int(end)}' for i, (start, end) in enumerate(range_list)])
    y_true_list = np.array(y_true_list)
    y_pred_list = np.array(y_pred_list)

    data = {
    'Batch': batch_list,
    'Y_True': y_true_list,
    'Y_Pred': y_pred_list,
    'Chain': chain_list,
    'Range': x_strings
    }

    results_df = pd.DataFrame(data)
    results_csv_path = f'{merged_config["result_dir"]}/{merged_config["prediction_csv"]}'
    results_df.to_csv(results_csv_path, index=False)
    print('results saved to csv:', results_csv_path)

if __name__ == "__main__":
    cluster = 'cluster1_8A_manual_rescale'
    spike_protein = 'Delta'
    proteins = [f'fold_vh16_vl106_seed2_model_0',
                f'fold_vh16_vl106_seed4_model_0',
                f'fold_vh16_vl106_seed4_model_1',
                f'fold_vh16_vl106_seed4_model_2',
                f'fold_vh16_vl106_seed5_model_2']
    proteins = [f'{spike_protein}_{protein}' for protein in proteins]
    graph_type = 'GearNet'
    epoch = [60, 70, 80, 90, 100]
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

    for protein in proteins:
        data_list = os.listdir(f'/home/lwang/models/HDX_LSTM/data/COVID_SPIKE/graph_ensemble_GearNet/{protein}/{cluster}')
        protein_list = [data.split('.')[0] for data in data_list]

        save_args = {
            #save setting
            'plot_title': 'test_set',
            'date': '240730',
            'model_name': 'GearNet',
            'version': 'v1',
            'result_dir': f'/home/lwang/models/HDX_LSTM/data/COVID_SPIKE/prediction/',
            'prediction_csv': None,
            'cluster':cluster,
            'protein_name': protein,
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
            # 240619 GearNet success in 1UGH, 6DJL, 6N1Z, 8F7A
            #save_args['model_path'] = f'/home/lwang/models/HDX_LSTM/results/240619_GearNetEdge/model_GVP_epoch60_cluster1_hop1_v{i}.pth'

            # 240730 GVP
            #save_args['model_path'] = f'/home/lwang/models/HDX_LSTM/results/240817_GVP/model_GVP56_cluster1_8A_manual_rad8+knn+seqEdge_layer5_v0.pth'
            #save_args['prediction_csv'] = f"HDX_pred_GVP56_{protein_name}_{cluster}_v{i}.csv"
            #main(save_args)

            #240920 GearNet better performance in test set
            '''save_args['model_path'] = f'/home/lwang/models/HDX_LSTM/results/240918_GVP/model_GN56_cluster1_8A_manual_rescale_v0_epoch{epoch[i]-1}.pth'
            save_args['prediction_csv'] = f"HDX_pred_GearNet56_{protein_name}_{cluster}_v{i}_repeat.csv"
            save_args['result_dir'] = f'/home/lwang/models/HDX_LSTM/data/hdock/prediction/GN56_epoch{epoch[i]}'
            if not os.path.exists(save_args['result_dir']):
                os.mkdir(save_args['result_dir'])
            main(save_args)'''

            #241128 GearNet-hse in covid spike protein prediction
            save_args['model_path'] = f'/home/lwang/models/HDX_LSTM/results/241110_GearNet/model_GN56_cluster1_8A_manual_rescale_v0_hse_epoch{epoch[i]-1}.pth'
            save_args['prediction_csv'] = f"HDX_pred_GearNet56_{protein}_epoch{epoch[i]}.csv"
            save_args['result_dir'] = f'/home/lwang/models/HDX_LSTM/data/COVID_SPIKE/predictions/{protein}'
            if not os.path.exists(save_args['result_dir']):
                os.makedirs(save_args['result_dir'])
            main(save_args)
    