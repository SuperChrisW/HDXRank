### training ###
import os
import sys
import copy

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.loader import DataLoader
from pepGraph_model import MixGCNBiLSTM, GCN
from pepGraph_BiLSTM import BiLSTM

import pandas as pd
import numpy as np
from tqdm import tqdm
from pdb2sql import interface

from scipy.stats import pearsonr
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
from pepGraph_utlis import seq_embedding, neighbor_search, neighbor_filter
from accelerate import Accelerator

def prepare_data(root_dir, embedding_file, pepGraph_fname, pdb, chain,
                 database_id, protein, state, correction_value=0):
    embedding_dir = os.path.join(root_dir, 'embedding_files')
    pepGraph_dir = os.path.join(root_dir, 'graph_ensemble')
    HDX_dir = os.path.join(root_dir, 'HDX_files')
    pdb_dir = os.path.join(root_dir, 'structure')

    embedding_file = f'{embedding_dir}/{embedding_file}'
    target_HDX_file = f'{HDX_dir}/{database_id}.xlsx'
    pepGraph_file = f'{pepGraph_dir}/{pepGraph_fname}'
    pdb_file = f'{pdb_dir}/{pdb}.pdb'

    if os.path.isfile(embedding_file) == False \
        or os.path.isfile(target_HDX_file) == False \
        or os.path.isfile(pepGraph_file) == False:
        #print('file not found')
        return False

    # sequence embedding process
    '''
    input_array, truth_array, start_pos, end_pos, log_t = \
                                        seq_embedding(target_HDX_file, embedding_file, pdb_file, protein, state, chain,
                                                        mismatch_report=False, correction_value=correction_value, filter = ['unique'])
    '''
    # graph process
    pepGraph_ensemble, index_dict = torch.load(pepGraph_file)

    return pepGraph_ensemble, index_dict
    #return input_array, truth_array, start_pos, end_pos, log_t, pepGraph

def data(hdx_df, root_dir):
    input_data = []
    progress_bar = tqdm(hdx_df.iterrows(), total=len(hdx_df))
    for row_id, row in hdx_df.iterrows():
        pdb = row['apo_identifier'].strip().split('.')[0]
        chain = row['chain_identifier']
        protein = row['protein']
        state = row['state']
        database_id = row['database_id']
        correction_value = row['correction_value']
        embedding_fname = f'{pdb}_{chain}.embedding.txt'
        pepGraph_fname = f'{pdb}.pt'
 
        print(database_id, pdb, chain)
        returned_value = prepare_data(root_dir, embedding_fname, pepGraph_fname, pdb, chain,
                                      database_id, protein, state, correction_value)
        
        if not returned_value:
            continue
        input_array, truth_array, start_pos, end_pos, log_t, pepGraph = returned_value

        protein_embedding = pepGraph['embedding']
        if np.isnan(pepGraph.x).any():
            print("pepGraph.x contains NaN values")
        elif len(start_pos) == 0:
            continue
        elif pepGraph.x.shape[0] != input_array.shape[0]:
            print('error: pepGraph.x and input_array have different length!!!!')
            continue

        pepGraph_data = []
        index_dict = pepGraph['pep_index'] ### format: dict[f'{chain}_{start_pos}_{end_pos}'] = {i}
        for seq_id, compressed_data in enumerate(zip(start_pos, end_pos, log_t)):
            start, end, timelog = compressed_data
            pep_label = f'{chain}_{start}_{end}'
            if pep_label in index_dict:
                graph = copy.deepcopy(pepGraph)
                graph['pep_index'] = int(index_dict[pep_label])
                graph.y = torch.tensor(truth_array[seq_id], dtype=torch.float32)
                
                # generate local graph consisted of two shells of neighbors
                neighbor1 = neighbor_search([graph['pep_index']], graph.edge_index.t())
                neighbors = set(list(neighbor1) + [graph['pep_index']])
                graph = neighbor_filter(graph, neighbors)
                
                #generate sequence embedding
                graph['embedding'] = input_array[graph['pep_index']]
                #graph['debug_label'] = f'{pdb}_{chain}_{start}_{end}'

                pepGraph_data.append(graph)
            else:
                print('pep_label not found:', pep_label)
                continue
        input_data.extend(pepGraph_data)
    return input_data

def train_model(model, num_epochs, optimizer, train_loader, val_loader, loss_fn):
    rp_train = []
    rmse_train_list = []
    rp_val = []
    rmse_val_list = []

    for epoch in range(num_epochs):
        list1_train = []
        list2_train = []
        epoch_train_losses = []

        list1_val = []
        list2_val = []
        epoch_val_losses = []
        model.train()
        for graph_batch in train_loader:
            graph_batch.to(device)
            targets = graph_batch.y.to(dtype=torch.float32)
            outputs = model(graph_batch)

            outputs = outputs.squeeze(-1).to(dtype = torch.float32)
            train_loss = loss_fn(outputs, targets)

            optimizer.zero_grad()
            train_loss.backward()
            optimizer.step()

            epoch_train_losses.append(train_loss.data.cpu().numpy())

            targets = targets.data.cpu().numpy()
            outputs = outputs.data.cpu().numpy()
            list1_train=np.append(list1_train,targets)
            list2_train=np.append(list2_train,outputs)

        model.eval()
        with torch.no_grad():
            val_loss = 0
            for graph_batch in val_loader:
                graph_batch.to(device)
                targets = graph_batch.y.to(dtype=torch.float32)
                outputs = model(graph_batch)

                outputs = outputs.squeeze(-1).to(dtype = torch.float32)
                outputs = outputs.to(dtype = torch.float32)
                val_loss = loss_fn(outputs, targets)

                targets=targets.data.cpu().numpy()
                outputs=outputs.data.cpu().numpy()
                list1_val=np.append(list1_val,targets)
                list2_val=np.append(list2_val,outputs)

                epoch_val_losses.append(val_loss.data.cpu().numpy())

        epoch_train_losses = np.mean(np.array(epoch_train_losses))
        epoch_val_losses = np.mean(np.array(epoch_val_losses))

        epoch_rp_train = np.corrcoef(list2_train, list1_train)[0,1]
        epoch_rp_val = np.corrcoef(list2_val, list1_val)[0,1]

        rp_train.append(np.mean(epoch_rp_train))
        rp_val.append(np.mean(epoch_rp_val))

        x = np.array(list1_train).reshape(-1,1)
        y = np.array(list2_train).reshape(-1,1)
        epoch_train_rmse= np.sqrt(((y - x) ** 2).mean())
        rmse_train_list.append(epoch_train_rmse)
        x = np.array(list1_val).reshape(-1,1)
        y = np.array(list2_val).reshape(-1,1)
        epoch_val_rmse= np.sqrt(((y - x) ** 2).mean())
        rmse_val_list.append(epoch_val_rmse)

        print('Epoch  Train Loss  Val Loss  PCC Train  PCC Val  Train RMSE  Val RMSE')
        print("{:5d}  {:10.3f}  {:9.3f}  {:9.3f}  {:8.3f}  {:10.3f}  {:9.3f}".format(
        epoch, epoch_train_losses, epoch_val_losses, epoch_rp_train, epoch_rp_val, epoch_train_rmse, epoch_val_rmse))

    return rmse_train_list, rmse_val_list, rp_train, rp_val

def test_model(model, test_loader, device):
    y_pred = []
    y_true = []
    model.eval()
    for i, graph_batch in enumerate(test_loader):
        graph_batch.to(device)
        targets = graph_batch.y.to(dtype=torch.float32)
        outputs = model(graph_batch)
        outputs = outputs.squeeze(-1).to(dtype = torch.float32)

        y_pred.append(outputs.cpu().detach().numpy())
        y_true.append(targets.cpu().detach().numpy())
    y_pred = np.concatenate(y_pred, axis=0)
    y_true = np.concatenate(y_true, axis=0)
    return y_true, y_pred

if __name__ == "__main__":

    ##################################### initial setting ##################################### 
    root_dir = "/home/lwang/AI-HDX-main/HDX_MS_dataset/database_collection/feature"
    summary_HDX_file = f'{root_dir}/merged_data.xlsx'
    hdx_df = pd.read_excel(summary_HDX_file, sheet_name='Sheet1')
    hdx_df = hdx_df.dropna(subset=['chain_identifier'])

    pepGraph_dir = os.path.join(root_dir, 'graph_ensemble')
    result_dir = '/home/lwang/models/HDX_LSTM/results/temp'
    if not os.path.exists(result_dir):
        os.mkdir(result_dir)
    ##################################### initial setting ##################################### 
    
    device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')

    ### data preparation ###
    apo_input = []
    complex_input = []
    hdx_df = hdx_df.drop_duplicates(subset=['apo_identifier'])
    for row_id, row in hdx_df.iterrows():
        pdb = row['apo_identifier'].strip().split('.')[0]
        pepGraph_file = f'{pepGraph_dir}/{pdb}.pt'
        if os.path.isfile(pepGraph_file):
            pepGraph_ensemble, index_dict = torch.load(pepGraph_file)
            if row['complex_state'] == 'single':
                apo_input.extend(pepGraph_ensemble)
            elif row['complex_state'] == 'complex':
                complex_input.extend(pepGraph_ensemble)
        else:
            continue

    ### model initialization ###
    args = {'in_dim': 48, 'final_hidden_dim': 48,
        'seq_dim': 10, 'struc_dim': 10, 'evo_dim': 10, 'time_dim': 5,
        'num_hidden_channels': 10, 'num_out_channels': 20, 
        'GNN_out_dim': 16, 'GNN_hidden_dim': 16,
        'drop_out': 0.3, 'num_GNN_layers': 2, 'GNN_type': 'GraphSAGE'}
    Mix_model = GCN(args).to(device)
    model = Mix_model

    ### training ###
    num_epochs = 1
    batch_size = 16
    num_workers = 4
    loss_fn = nn.BCELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001, weight_decay=5e-4)

    train_apo, val_apo = train_test_split(apo_input, test_size=0.3, random_state=42)
    train_complex, val_complex = train_test_split(complex_input, test_size=0.3, random_state=42)

    train_set = train_apo + train_complex
    val_set = val_apo + val_complex
    train_loader = DataLoader(train_set, batch_size = batch_size, shuffle=True, num_workers=num_workers)
    val_loader =  DataLoader(val_set, batch_size = batch_size, shuffle=False, num_workers=num_workers)
    print('length of train_Set:', len(train_set))
    print('length of val_Set:', len(val_set))

    rmse_train_list, rmse_val_list, rp_train, rp_val = train_model(
        model, num_epochs, optimizer, train_loader, val_loader, loss_fn)

    ### save model ###
    torch.save(model, f'{result_dir}/model.pt')

    # PCC Plot
    plt.figure(figsize=(10, 6))
    plt.plot(rp_train, label='Training PCC')
    plt.plot(rp_val, label = 'Val PCC')
    plt.title('Pearson Correlation Coefficient over epochs')
    plt.ylim(0, 1)
    plt.xlabel('Epochs')
    plt.ylabel('PCC')
    plt.legend()
    plt.savefig(f'{result_dir}/pcc_plot.png')  # Save the plot as a PNG file

    # RMSE Plot
    plt.figure(figsize=(10, 6))
    plt.plot(rmse_train_list, label='train RMSE')
    plt.plot(rmse_val_list, label='val RMSE')
    plt.title('Root Mean Square Error over epochs')
    plt.ylim(0.05, 0.2)
    plt.xlabel('Epochs')
    plt.ylabel('RMSE')
    plt.legend()
    plt.savefig(f'{result_dir}/rmse_plot.png')  # Save the plot as a PNG file

    y_true, y_pred = test_model(model, val_loader, device)
    plt.figure()
    plt.scatter(y_true, y_pred, alpha=0.5)
    plt.xlabel('Experimental HDX')
    plt.ylabel('Predicted HDX')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.plot([0, 1], [0, 1], '--', color='gray')
    plt.savefig(f'{result_dir}/val_result.png')
