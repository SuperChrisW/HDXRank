### prediction ###
import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.loader import DataLoader
from pepGraph_model import MixGCNBiLSTM
from pepGraph_BiLSTM import BiLSTM

import pandas as pd
import numpy as np

from scipy.stats import pearsonr
import matplotlib.pyplot as plt

from pepGraph_train import data
from pepGraph_basemodel_training import base_data, test_basemodel

def test_pepGraph(model, test_loader, device):
    y_pred = []
    y_true = []
    pep_centers = []
    model.eval()
    for i, graph_batch in enumerate(test_loader):
        graph_batch.to(device)
        embedding_input = graph_batch['embedding']
        targets = graph_batch.y.to(dtype=torch.float64)
        embedding_input = embedding_input.unsqueeze(1)
        outputs = model(graph_batch, embedding_input)
        outputs = outputs.squeeze(-1).to(dtype = torch.float64)

        positions = graph_batch['range']
        for pos in positions:
            pep_centers.append((pos[1].item() + pos[0].item()) / 2)

        y_pred.append(outputs.cpu().detach().numpy())
        y_true.append(targets.cpu().detach().numpy())
    y_pred = np.concatenate(y_pred, axis=0)
    y_true = np.concatenate(y_true, axis=0)
    pep_centers = torch.tensor(pep_centers, dtype=torch.float).cpu().detach().numpy()
    return y_true, y_pred, pep_centers

if __name__ == "__main__":
    ##################################### initial setting #####################################
    root_dir = "/home/lwang/AI-HDX-main/HDX_MS_dataset/complexpair_dataset"
    HDX_summary_file = f'{root_dir}/merged_complex_pair.xlsx'
    hdx_df = pd.read_excel(HDX_summary_file, sheet_name='Sheet1')
    hdx_df = hdx_df.dropna(subset=['chain_identifier'])

    dssp_dir = os.path.join(root_dir, 'dssp_files')
    hhm_dir = os.path.join(root_dir, 'hhm_files')
    rigidity_dir = os.path.join(root_dir, 'rigidity_files')
    embedding_dir = os.path.join(root_dir, 'embedding_files')
    proflex_dir = os.path.join(root_dir, 'proflex_files')
    pepGraph_dir = os.path.join(root_dir, 'pepGraph_files')

    base_model_path = '/home/lwang/models/HDX_LSTM/results/temp_basemodel/model.pt'
    model_path = '/home/lwang/models/HDX_LSTM/results/240228_pepGraph_newfeature/model.pt'

    ##################################### initial setting #####################################

    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    print(device)    

    choice = 'pepGraph'
    save_dir = '/home/lwang/models/HDX_LSTM/results/240228_pepGraph_newfeature'
    ### data preparation ###    
    if choice == 'pepGraph':
        input_data = data(hdx_df, root_dir)
        print('length of input_data:', len(input_data))

        Mix_model = torch.load(model_path, map_location=device)
        model = Mix_model
        model.to(device)

        batch_size = 64
        test_loader = DataLoader(input_data, batch_size = batch_size, shuffle=False)
        y_true, y_pred, pep_centers = test_pepGraph(model, test_loader, device)

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
        y_true, y_pred = test_basemodel(model, test_loader, device)
        
    pcc = pearsonr(y_true, y_pred)
    print('PCC:', pcc[0])

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
    plt.legend(['PCC: {:.2f}'.format(pcc[0])])
    plt.savefig(f'{save_dir}/temp_test_contact484.png')