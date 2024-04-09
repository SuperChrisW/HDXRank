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

from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt

#from pepGraph_train import data
#from pepGraph_basemodel_training import base_data, test_basemodel
from sklearn.metrics import mean_absolute_error, mean_squared_error
from math import sqrt

def test_model(model, test_loader, device):
    y_pred = []
    y_true = []
    range_list = []
    model.eval()
    for i, graph_batch in enumerate(test_loader):
        graph_batch.to(device)
        targets = graph_batch.y.to(dtype=torch.float32)
        outputs = model(graph_batch)
        outputs = outputs.to(dtype = torch.float32)
        range_list.extend(graph_batch.range.cpu().detach().numpy())
        y_pred.append(outputs.cpu().detach().numpy())
        y_true.append(targets.cpu().detach().numpy())
    y_pred = np.concatenate(y_pred, axis=0)
    y_true = np.concatenate(y_true, axis=0)
    return y_true, y_pred, range_list

if __name__ == "__main__":
    ##################################### initial setting #####################################
    root_dir = "/home/lwang/models/HDX_LSTM/data/COVID_SPIKE"
    HDX_summary_file = f'{root_dir}/COVID_record.xlsx'
    hdx_df = pd.read_excel(HDX_summary_file, sheet_name='Sheet1')
    hdx_df = hdx_df.dropna(subset=['chain_identifier'])

    dssp_dir = os.path.join(root_dir, 'dssp_files')
    hhm_dir = os.path.join(root_dir, 'hhm_files')
    rigidity_dir = os.path.join(root_dir, 'rigidity_files')
    embedding_dir = os.path.join(root_dir, 'embedding_files')
    proflex_dir = os.path.join(root_dir, 'proflex_files')
    pepGraph_dir = os.path.join(root_dir, 'graph_ensemble', 'v2_ensemble')

    model_path = f'/home/lwang/models/HDX_LSTM/results/240402_GAT_v2/hop1/best_model.pt'
    ##################################### initial setting #####################################

    accelerator = Accelerator()
    device = accelerator.device
    #device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

    choice = 'pepGraph'
    ### data preparation ###    

    input_data = []
    hdx_df = hdx_df.drop_duplicates(subset=['apo_identifier'])
    for row_id, row in hdx_df.iterrows():
        pdb = row['apo_identifier'].strip().split('.')[0]
        if not pdb == 'Omicron_spike':
            continue
        pepGraph_file = f'{pepGraph_dir}/{pdb}.pt'
        if os.path.isfile(pepGraph_file):
            pepGraph_ensemble = torch.load(pepGraph_file) # list of graphs
            input_data.extend(pepGraph_ensemble)
        else:
            continue

    print('length of input_data:', len(input_data))
    if choice == 'pepGraph':
        config = {
            'num_epochs': 120,
            'batch_size': 16,
            'learning_rate': 0.001,
            'weight_decay': 5e-4,
            'GNN_type': 'GAT',
            'num_GNN_layers': 3,
            'cross_validation_num': 1,
            'num_workers': 4
        }
        training_args = {'in_dim': 91, 'final_hidden_dim': 48,
                'seq_dim': 10, 'struc_dim': 10, 'evo_dim': 10, 'time_dim': 5,
                'num_hidden_channels': 10, 'num_out_channels': 20, 

                'feat_in_dim': 45, 'topo_in_dim': 42, 'num_heads': 8,
                'GNN_out_dim': 1, 'GNN_hidden_dim': 32,
                'drop_out': 0.5, 'num_GNN_layers': config['num_GNN_layers'], 'GNN_type': config['GNN_type'],
                'graph_hop': 'hop1',
                'result_dir': '/home/lwang/models/HDX_LSTM/data/COVID_SPIKE'
        }
        model = GNN_v2(training_args)
        unwrapped_model = accelerator.unwrap_model(model)
        path_to_checkpoint = os.path.join(model_path, "pytorch_model.bin")
        unwrapped_model.load_state_dict(torch.load(path_to_checkpoint))
        unwrapped_model.to(device)

        batch_size = config['batch_size']
        test_loader = DataLoader(input_data, batch_size = batch_size, shuffle=False)
        test_loader = accelerator.prepare(test_loader)
        y_true, y_pred, range_list = test_model(model, test_loader, device)

    '''
    ### evaluation ###
    y_pred = y_pred.flatten()
    print(y_true.shape, y_pred.shape)
    pcc = pearsonr(y_true, y_pred)
    mae = mean_absolute_error(y_true, y_pred)
    rmse = sqrt(mean_squared_error(y_true, y_pred))
    spR = spearmanr(y_true, y_pred)
    print('Mean Absolute Error (MAE):', mae)
    print('PCC:', pcc[0])    
    print('Root Mean Squared Error (RMSE):', rmse)
    print('Spearman:', spR[0])

    ### plot ###
    plt.figure()
    plt.scatter(y_true, y_pred, alpha=0.5)
    plt.xlabel('Experimental HDX')
    plt.ylabel('Predicted HDX')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.plot([0, 1], [0, 1], '--', color='gray')  # Add diagonal dashed line
    plt.scatter(y_true, y_pred, alpha=0.5)
    plt.xlabel('Experimental HDX')
    plt.ylabel('Predicted HDX')
    plt.legend(['PCC: {:.3f}'.format(pcc[0]), 'MAE: {:.3f}'.format(mae), 'RMSE: {:.3f}'.format(rmse), 'Spearman: {:.3f}'.format(spR[0])])

    date = '240408'
    model_name = 'GAT'
    version = 'v2.1'
    dataset = 'COVID_SPIKE'
    save_dir = training_args['result_dir']
    plt.savefig(f'{save_dir}/{date}_{model_name}_{version}_{dataset}.png')
    '''
    range_list = np.array(range_list).reshape(-1, 2)
    x_label = np.mean(range_list, axis=1)
    print(x_label.shape, y_true.shape, y_pred.shape)

    def slide_mean(data, window_size=5):
        return np.convolve(data, np.ones(window_size)/window_size, mode='same')
    x_label = slide_mean(x_label, window_size=3)
    y_true = slide_mean(y_true.flatten(), window_size=3)
    y_pred = slide_mean(y_pred.flatten(), window_size=3)

    x_label, y_true, y_pred = zip(*sorted(zip(x_label, y_true, y_pred)))

    plt.figure(figsize=(10, 6))  # Adjust figure size for better visibility
    plt.plot(x_label, y_true, label='True HDX', color='green', marker = 'o', linestyle='-', linewidth=2, markersize=3)
    plt.plot(x_label, y_pred, label='Predicted HDX', color='blue', marker='s', linestyle='--', linewidth=2, markersize=3)

    plt.title('Comparison of True and Predicted HDX: Omicron spike protein')  # Adding chart title
    plt.xlabel('Peptide Center (Residue Index)')  # Labeling the x-axis
    plt.ylabel('HDX Level')  # Labeling the y-axis
    plt.legend(loc='best')  # Displaying the legend in the best location
    plt.grid(True)  # Adding gridlines for better readability
    plt.savefig(f'/home/lwang/models/HDX_LSTM/data/COVID_SPIKE/Omicronspike_linechart2_.png')  # Save the plot as a PNG file