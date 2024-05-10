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
from pepGraph_BiLSTM import BiLSTM

import pandas as pd
import numpy as np

from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt

from pepGraph_train import data
from pepGraph_basemodel_training import base_data, test_basemodel
from sklearn.metrics import mean_absolute_error, mean_squared_error
from math import sqrt

def test_model(model, test_loader, device):
    y_pred = []
    y_true = []
    range_label = []
    model.eval()
    for i, graph_batch in enumerate(test_loader):
        graph_batch.to(device)
        targets = graph_batch.y.to(dtype=torch.float32)
        outputs = model(graph_batch)
        outputs = outputs.to(dtype = torch.float32)

        range_label.extend(graph_batch.label)
        y_pred.append(outputs.cpu().detach().numpy())
        y_true.append(targets.cpu().detach().numpy())
    y_pred = np.concatenate(y_pred, axis=0)
    y_true = np.concatenate(y_true, axis=0)

    return y_true, y_pred, range_label

def predict(input_data, model_path, accelerator, choice = 'pepGraph'):
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
        args = {'in_dim': 48, 'final_hidden_dim': 48,
                'seq_dim': 10, 'struc_dim': 10, 'evo_dim': 10, 'time_dim': 5,
                'num_hidden_channels': 10, 'num_out_channels': 20, 
                'GNN_out_dim': 16, 'GNN_hidden_dim': 16,
                'drop_out': 0.5, 'num_GNN_layers': config['num_GNN_layers'], 'GNN_type': config['GNN_type']
        }
        model = MixGCNBiLSTM(args)
        device = accelerator.device

        unwrapped_model = accelerator.unwrap_model(model)
        path_to_checkpoint = os.path.join(model_path, "pytorch_model.bin")
        unwrapped_model.load_state_dict(torch.load(path_to_checkpoint))
        #model.gcn.GNN_type = 'GAT'
        #model.GNN_type = 'GAT'
        unwrapped_model.to(device)

        batch_size = 16
        test_loader = DataLoader(input_data, batch_size = batch_size, shuffle=False)
        test_loader = accelerator.prepare(test_loader)
        y_true, y_pred, range_label = test_model(unwrapped_model, test_loader, device)

    elif choice == 'base_model':
        input_data, truth_data = base_data(hdx_df, root_dir)
        input_data = data(hdx_df, f'/home/lwang/AI-HDX-main/HDX_MS_dataset/database_collection/feature')
        print('length of input_data:', len(input_data))
        input_data = torch.tensor(input_data, dtype=torch.float32)
        truth_data = torch.tensor(truth_data, dtype=torch.float32)

        base_model = torch.load(base_model_path, map_location=device)
        model = base_model
        model.to(device)

        batch_size = 64
        test_loader = DataLoader(list(zip(input_data, truth_data)), batch_size = batch_size, shuffle=False)
        test_loader = accelerator.prepare(test_loader)
        y_true, y_pred = test_basemodel(model, test_loader, device)

    return y_true, y_pred, range_label

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
    pepGraph_dir = os.path.join(root_dir, 'graph_ensemble/v1_ensemble')

    base_model_path = '/home/lwang/models/HDX_LSTM/results/temp_basemodel/model.pt'
    model_path = f'/home/lwang/models/HDX_LSTM/results/240327_BiLSTM_GAT/best_model.pt'
    save_dir = '/home/lwang/models/HDX_LSTM/data/COVID_SPIKE'
    accelerator = Accelerator()    
    #device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')    
    ##################################### initial setting #####################################
    ### data preparation ### 
    target_protein = 'Wuhan_spike'
    target_chain = ['A']
    Truth_list = []
    Pred_list = []
    Range_list = []

    for chain in target_chain:
        input_data = []
        #hdx_df = hdx_df.drop_duplicates(subset=['apo_identifier'])
        temp_df = hdx_df[(hdx_df['apo_identifier'] == f'{target_protein}.pdb') & (hdx_df['chain_identifier'] == chain)]

        for row_id, row in temp_df.iterrows():
            pdb = row['apo_identifier'].strip().split('.')[0]
            pepGraph_file = f'{pepGraph_dir}/{pdb}.pt'
            if os.path.isfile(pepGraph_file):
                pepGraph_ensemble, _ = torch.load(pepGraph_file)
                input_data.extend(pepGraph_ensemble)
            else:
                continue
        print('length of input_data:', len(input_data))
        y_true, y_pred, range_label = predict(input_data, model_path, accelerator)
        Truth_list.append(y_true)
        Pred_list.append(y_pred)
        Range_list.append(range_label)

    ### evaluation: line chart residue graph ###

    plt.figure(figsize=(10, 6))
    colors = ['blue', 'green', 'red']
    markers = ['o', 's', '^']

    def slide_mean(data, window_size=3):
        return np.convolve(data, np.ones(window_size)/window_size, mode='same')
    for i, (y_true, y_pred, range_label) in enumerate(zip(Truth_list, Pred_list, Range_list)):
        range_list = []
        for label in range_label:
            parts = label.split('_')
            range_list.append(int(parts[1]))
            range_list.append(int(parts[2]))
        range_list = np.array(range_list).reshape(-1, 2)
        x_label = np.mean(range_list, axis=1)
        y_true = y_true.flatten()
        y_pred = y_pred.flatten()

        x_label, y_true, y_pred = zip(*sorted(zip(x_label, y_true, y_pred)))

        y_true = slide_mean(y_true, window_size=3)
        y_pred = slide_mean(y_pred, window_size=3)
        x_label = slide_mean(x_label, window_size=3)
        x_label = x_label[1:-1]
        y_true = y_true[1:-1]
        y_pred = y_pred[1:-1]

        plt.plot(x_label, y_pred, label=f'Predicted Chain {target_chain[i]}', color=colors[i], marker=markers[i], linestyle='--', linewidth=2, markersize=3)

    plt.plot(x_label, y_true, label='True HDX', color='grey', marker='o', linestyle='-', linewidth=2, markersize=3)
    plt.axvspan(330, 530, color='yellow', alpha=0.2, label = 'Receptor binidng region (RBD)')

    plt.title('Comparison of True and Predicted HDX: Delta spike protein')  # Adding chart title
    plt.xlabel('Peptide Center (Residue Index)')  # Labeling the x-axis
    plt.ylabel('HDX Level')  # Labeling the y-axis
    plt.legend(loc='best')  # Displaying the legend in the best location
    plt.grid(True)  # Adding gridlines for better readability

    date = '240408'
    model_name = 'BiLSTM_GAT'
    version = 'v1'
    dataset = 'Delta_linechart_'
    plt.savefig(f'{save_dir}/{date}_{model_name}_{version}_{dataset}.png')


    '''       
    ### evaluation: scatter plot ###
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
    '''