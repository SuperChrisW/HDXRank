### prediction ###
import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

from accelerate import Accelerator
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.loader import DataLoader
from pepGraph_model import MixGCNBiLSTM
from pepGraph_BiLSTM import GNN_v2

import pandas as pd
import numpy as np
import itertools

from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt
from sklearn.metrics import mean_absolute_error, mean_squared_error
from math import sqrt

def test_model(model, test_loader, device):
    y_pred = []
    y_true = []
    range_list = []
    chain_list = []
    model.eval()
    for i, graph_batch in enumerate(test_loader):
        graph_batch.to(device)
        targets = graph_batch.y.to(dtype=torch.float32) #only for Full set and RTT-BCD set
        outputs = model(graph_batch)
        outputs = outputs.to(dtype = torch.float32)
        range_list.extend(graph_batch.range.cpu().detach().numpy())
        chain_list.extend(graph_batch.chain)

        y_pred.append(outputs.cpu().detach().numpy())
        y_true.append(targets.cpu().detach().numpy())
    y_pred = np.concatenate(y_pred, axis=0)
    y_true = np.concatenate(y_true, axis=0)
    return y_true, y_pred, range_list, chain_list

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

def run_prediction(model, hdx_df, merged_config, pepGraph_dir, device):
    ##################################### data preparation #####################################
    test_loader = load_data(hdx_df, pepGraph_dir, merged_config)
    print('data loaded successfully!')
    y_true, y_pred, range_list, chain_list = test_model(model, test_loader, device)
    return y_true, y_pred, range_list, chain_list

def load_model(model_path, config, accelerator):
    model = MixGCNBiLSTM(config)
    unwrapped_model = accelerator.unwrap_model(model)
    path_to_checkpoint = os.path.join(model_path, "pytorch_model.bin")
    unwrapped_model.load_state_dict(torch.load(path_to_checkpoint))
    return unwrapped_model

def load_data(df, pepGraph_dir, merged_config):
    input_data = []
    temp_df = df.drop_duplicates(subset=['apo_identifier'], keep='first')
    load_proteins = merged_config['load_proteins'] if isinstance(merged_config['load_proteins'], list) else [merged_config['load_proteins']]

    for row_id, row in temp_df.iterrows():
        pdb = row['apo_identifier'].strip().split('.')[0]
        if pdb in load_proteins:
            pepGraph_file = os.path.join(pepGraph_dir, f'{pdb}.pt')
            if os.path.isfile(pepGraph_file):
                pepGraph_ensemble = torch.load(pepGraph_file) # list of graphs
                input_data.extend(pepGraph_ensemble)
    print('length of input_data:', len(input_data))
    return DataLoader(input_data, batch_size=merged_config['batch_size'], shuffle=True)

def finetune(model, finetune_loader, config, accelerator):
    
    '''
    # freeze some layers
    for param in model.parameters():
        param.requires_grad = False
    model.fc1.requires_grad = True
    model.fc2.requires_grad = True 
    optimizer = torch.optim.Adam([
        {'params': model.fc1.parameters()},
        {'params': model.fc2.parameters()}
    ], lr=config['learning_rate'])
    '''

    loss_fn = torch.nn.BCELoss()  # Assuming a regression problem; change if necessary
    optimizer = torch.optim.Adam(model.parameters(), lr=config['learning_rate'], weight_decay=config['weight_decay'])

    model.train()
    for epoch in range(config['num_epochs']):
        for graph_batch in finetune_loader:
            targets = graph_batch.y.to(dtype=torch.float32)
            outputs = model(graph_batch)
            train_loss = loss_fn(outputs, targets)

            optimizer.zero_grad()
            accelerator.backward(train_loss)
            optimizer.step()
        print(f'Epoch {epoch+1}, Loss: {train_loss.item()}')
    return model

def main(save_args):
    ##################################### initial setting #####################################
    root_dir = "/home/lwang/models/HDX_LSTM/data/COVID_SPIKE"
    HDX_summary_file = f'{root_dir}/hdock/COVID_Delta.xlsx'
    hdx_df = pd.read_excel(HDX_summary_file, sheet_name='Sheet1')
    hdx_df = hdx_df.dropna(subset=['chain_identifier'])
    pepGraph_dir = os.path.join(root_dir, 'graph_ensemble')

    model_path = save_args['model_path']

    accelerator = Accelerator()
    device = accelerator.device
    #device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    ##################################### config setting #####################################
    config = {
    #training setting
    'num_epochs': 8, # epochs for finetuning
    'batch_size': 16,
    'learning_rate': 0.001,
    'weight_decay': 5e-4,
    'GNN_type': 'GAT',
    'num_GNN_layers': 3,
    'cross_validation_num': 1,
    'num_workers': 4,
    'result_dir': '/home/lwang/models/HDX_LSTM/results/240411_BiLSTMGAT_v2_1',
    'data_log': False,  

    #model setting
    'num_hidden_channels': 10,
    'num_out_channels': 20,
    'feat_in_dim': 44,
    'topo_in_dim': 42,
    'num_heads': 8,
    'GNN_hidden_dim': 32,
    'GNN_out_dim': 16,
    'LSTM_out_dim': 16,
    'final_hidden_dim': 16,
    'drop_out': 0.5,
    'graph_hop': 'hop1'
    }
    merged_config = {**config, **save_args}

    model = load_model(model_path, merged_config, accelerator)
    model = model.to(device)
    # finetune here
    if merged_config['finetune_check']:
        merged_config['load_proteins'] = merged_config['finetune_protein']
        finetune_loader = load_data(hdx_df, pepGraph_dir, merged_config)
        finetune_loader = accelerator.prepare(finetune_loader)
        model = finetune(model, finetune_loader, merged_config, accelerator)

    # run_prediciton here 
    '''# plot RMSE vs seq_lens
    merged_config['load_proteins'] = merged_config['prediction_protein']
    y_true, y_pred, range_, chain = run_prediction(model, hdx_df, merged_config, pepGraph_dir, device)

    range_ = np.array(range_).reshape(-1, 2)
    seq_lens = (range_[:, 1] - range_[:, 0]).flatten() + 1


    mse_perSeqlen = {}
    for len in np.unique(seq_lens):
        mask = (seq_lens == len)
        rmse = sqrt(mean_squared_error(y_true[mask], y_pred[mask]))
        var = np.var(y_true[mask] - y_pred[mask])
        mse_perSeqlen[len] = rmse
        print(f'Sequence Length {len} MSE:', rmse , 'Variance:', var)
    
    #plot a graph of mse vs seq_lens bar chart
    plt.figure(figsize=(10, 5))
    plt.bar(mse_perSeqlen.keys(), mse_perSeqlen.values(), color='b')
    plt.xlabel('Sequence Length')
    plt.ylabel('RMSE')
    plt.title('RMSE vs Sequence Length')
    plt.savefig(f'{merged_config["result_dir"]}/RMSE_vs_Sequence_Length.png')
    '''
    
    y_true_list = []
    y_pred_list = []
    range_list = []
    batch_list = []
    chain_list = []

    prediction_list = merged_config['prediction_protein'] + merged_config['finetune_protein']
    for index, protein in enumerate(prediction_list):
        merged_config['load_proteins'] = protein
        print(protein)
        y_true, y_pred, range_, chain = run_prediction(model, hdx_df, merged_config, pepGraph_dir, device)

        if protein == 'bcd_6NZ2':
            chain = ['B'] * len(y_true)
        if protein in merged_config['finetune_protein']:
            batch = np.full(len(y_true), f'finetune')
        else:
            batch = np.full(len(y_true), protein)
        batch_list.extend(batch)
        y_true_list = np.concatenate((y_true_list, y_true), axis=0)
        y_pred_list = np.concatenate((y_pred_list, y_pred), axis=0)
        chain_list = np.concatenate((chain_list, chain), axis=0)
        range_list = np.concatenate((range_list, range_), axis=0)

    #plot_result(y_true_list, y_pred_list, range_list, batch_list, chain_list, merged_config)

    range_list = np.array(range_list).reshape(-1, 2)
    chain_list = np.array(chain_list)
    x_strings = np.array([f'{chain}_{int(start)}-{int(end)}' for i, (start, end) in enumerate(range_list) for chain in chain_list[i]])

    print(len(batch_list))
    print(y_true_list.shape)
    print(y_pred_list.shape)
    print(x_strings.shape)
    print(chain_list.shape)
    data = {
    'Batch': batch_list,
    'Y_True': y_true_list,
    'Y_Pred': y_pred_list,
    'Chain': chain_list,
    'Range': x_strings #FIXME: change to start, end
        # FIXME: add 'exposure', 'maxuptake', 'uptake'
    }
    results_df = pd.DataFrame(data)
    results_csv_path = f'{merged_config["result_dir"]}/{merged_config["prediction_csv"]}'
    results_df.to_csv(results_csv_path, index=False)
    print('results saved to csv:', results_csv_path)

    '''
    batch_list = np.array(batch_list)
    range_list = np.array(range_list).reshape(-1, 2)
    pepcenter = np.mean(range_list, axis=1)
    label_list = np.array([f'{chain}_{int(start)}-{int(end)}' for i, (start, end) in enumerate(range_list) for chain in chain_list[i]])
    Y_diff = y_true_list - y_pred_list

    base_mask = (batch_list == 'finetune')
    Y_diff_baseline = Y_diff[base_mask]
    base_label = label_list[base_mask]
    base_center = pepcenter[base_mask]
    Y_diff_baseline, base_label, base_center = zip(*sorted(zip(Y_diff_baseline, base_label, base_center)))

    for batch in np.unique(batch_list):
        batch_mask = (batch_list == batch)
        Y_diff_batch = Y_diff[batch_mask]
        batch_label = label_list[batch_mask]

        label_mask = np.isin(batch_label, base_label)
        Y_diff_batch = Y_diff_batch[label_mask]
        batch_label = batch_label[label_mask]
        label_mask = np.isin(base_label, batch_label)
        temp_Y_diff_baseline = np.array(Y_diff_baseline)[label_mask]

        # align batch_label and base_label
        batch_center = pepcenter[batch_mask]
        batch_center, Y_diff_batch, batch_label = zip(*sorted(zip(batch_center, Y_diff_batch, batch_label)))

        diff_rmse = mean_absolute_error(Y_diff_batch, temp_Y_diff_baseline)
        print(f'{batch} MAE:', diff_rmse)
    '''

if __name__ == "__main__":
    protein_name  = 'Hdock_DeltaRBD+V16noext'

    protein_list = [f'{protein_name}_{i}_revised' for i in range(1, 101)]

    #data_list = os.listdir('/home/lwang/models/S3214dataset/graph_ensemble')
    #protein_list = [data.split('.')[0] for data in data_list]

    save_args = {
        #save setting
        'plot_title': 'Hdock_AF_rtt+bcd_G0',
        'date': '240422',
        'model_name': 'BiLSTMGAT_finetune',
        'version': 'v2',
        'result_dir': f'/home/lwang/models/HDX_LSTM/results/240411_BiLSTMGAT_v2_1',
        'prediction_csv': f'HDX_pred_{protein_name}.csv',

        #prediction setting
        'prediction_protein': protein_list,
        'finetune_check' : False,
        'finetune_protein': ['delta_spikeprotein'],
        'model_path': '/home/lwang/models/HDX_LSTM/results/240411_BiLSTMGAT_v2_1/best_model.pt',

        #plot setting
        'slide_window': 1,
        'show_truth': True,
        'show_pred': True,
        }
    
    main(save_args)
    