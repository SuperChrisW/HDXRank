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
from pepGraph_BiLSTM import BiLSTM

import pandas as pd
import numpy as np
from tqdm import tqdm

from scipy.stats import pearsonr
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
from pepGraph_utlis import seq_embedding
from pdb2sql import interface

def prepare_data(root_dir, embedding_file, pdb, chain,
                 database_id, protein, state, correction_value=0):
    
    embedding_dir = os.path.join(root_dir, 'embedding_files')
    HDX_dir = os.path.join(root_dir, 'HDX_files')
    pdb_dir = os.path.join(root_dir, 'structure')

    embedding_file = f'{embedding_dir}/{embedding_file}'
    #target_HDX_file = f'{HDX_dir}/{database_id}.xlsx'
    target_HDX_file = f'{HDX_dir}/COVID_SPIKE.xlsx'
    pdb_file = f'{pdb_dir}/{pdb}.pdb'

    if os.path.isfile(embedding_file) == False \
        or os.path.isfile(target_HDX_file) == False:
        print('file not found')
        return False
    
    #db = interface(pdb_file) # pdb2sql interface object 
    db = None
    input_array, truth_array, start_pos, end_pos, log_t = \
                                        seq_embedding(target_HDX_file, embedding_file, pdb_file, protein, state, chain,
                                                        mismatch_report=True, correction_value=correction_value, filter = ['unique'])
    range_list = np.column_stack([start_pos, end_pos])
    return input_array, truth_array, range_list

def base_data(hdx_df, root_dir):
    input_apo = []
    truth_apo = []
    input_complex = []
    truth_complex = []
    range_list = []
    for row_id, row in hdx_df.iterrows():
        pdb = row['apo_identifier'].strip().split('.')[0]
        chain = row['chain_identifier']
        protein = row['protein']
        state = row['state']
        database_id = row['database_id']
        correction_value = row['correction_value']
        embedding_fname = f'{pdb}_{chain}.embedding.txt'
        complex_state = row['complex_state']
 
        print(database_id, pdb, chain)
        returned_value = prepare_data(root_dir, embedding_fname, pdb, chain,
                                      database_id, protein, state, correction_value)
        if not returned_value:
            continue
        input_array, truth_array, pos = returned_value
        range_list.extend(pos)
        
        if np.isnan(truth_array).any():
            print('nan value found')
            continue
        truth_array = np.array(truth_array)/100
        
        if complex_state == 'single':
            input_apo.extend(input_array)
            truth_apo.extend(truth_array)
        elif complex_state == 'complex':
            input_complex.extend(input_array)
            truth_complex.extend(truth_array)

    input_apo = np.array(input_apo)
    truth_apo = np.array(truth_apo)
    input_complex = np.array(input_complex)
    truth_complex = np.array(truth_complex)
    print('input_apo:', input_apo.shape)
    print('truth_apo:', truth_apo.shape)
    print('input_complex:', input_complex.shape)
    print('truth_complex:', truth_complex.shape)
    return input_apo, truth_apo, input_complex, truth_complex, range_list

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
        for data in train_loader:
            inputs, targets = data
            inputs, targets = inputs.to(device), targets.to(device)

            inputs = inputs.unsqueeze(1)
            outputs = model(inputs)

            outputs = outputs.squeeze(-1)
            outputs = outputs.to(dtype = torch.float32)
            train_loss = loss_fn(outputs, targets)

            optimizer.zero_grad()  # Clear existing gradients
            train_loss.backward()  # Compute gradients
            optimizer.step()  # Update weights

            epoch_train_losses.append(train_loss.data.cpu().numpy())

            targets = targets.data.cpu().numpy()
            outputs = outputs.data.cpu().numpy()
            list1_train=np.append(list1_train,targets)
            list2_train=np.append(list2_train,outputs)

        model.eval()
        with torch.no_grad():
            val_loss = 0
            for data in val_loader:
                inputs, targets = data
                inputs, targets = inputs.to(device), targets.to(device)
                
                inputs = inputs.unsqueeze(1)
                outputs = model(inputs)

                outputs = outputs.squeeze(-1)
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
        #print(list2_val, list1_val)

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

def test_basemodel(model, test_loader, device):
    y_pred = []
    y_true = []
    pep_centers = []
    model.eval()
    for i, data in enumerate(test_loader):
        inputs, targets = data
        inputs, targets = inputs.to(device), targets.to(device)
        
        inputs = inputs.unsqueeze(1)
        outputs = model(inputs)
        outputs = outputs.squeeze(-1)

        y_pred.append(outputs.cpu().detach().numpy())
        y_true.append(targets.cpu().detach().numpy())
    y_pred = np.concatenate(y_pred, axis=0)
    y_true = np.concatenate(y_true, axis=0)
    return y_true, y_pred

if __name__ == "__main__":
    ##################################### initial setting ##################################### 
    root_dir = "/home/lwang/AI-HDX-main/HDX_MS_dataset/database_collection/feature"
    summary_HDX_file = f'{root_dir}/merged_apo.xlsx'
    hdx_df = pd.read_excel(summary_HDX_file, sheet_name='Sheet1')
    hdx_df = hdx_df.dropna(subset=['chain_identifier'])

    dssp_dir = os.path.join(root_dir, 'dssp_files')
    hhm_dir = os.path.join(root_dir, 'hhm_files')
    rigidity_dir = os.path.join(root_dir, 'rigidity_files')
    proflex_dir = os.path.join(root_dir, 'proflex_files')
    embedding_dir = os.path.join(root_dir, 'embedding_files')
    HDX_dir = os.path.join(root_dir, 'HDX_files')

    result_dir = '/home/lwang/models/HDX_LSTM/results/temp_basemodel'
    if not os.path.exists(result_dir):
        os.mkdir(result_dir)
    ##################################### initial setting ##################################### 
        
    device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')
    print(device)

    ### data preparation ###
    input_data, truth_data = base_data(hdx_df, root_dir)
    print('length of input_data:', len(input_data))
    input_data = torch.tensor(input_data, dtype=torch.float32)
    truth_data = torch.tensor(truth_data, dtype=torch.float32)

    ### model initialization ###
    BiLSTM_model = BiLSTM().to(device)
    model = BiLSTM_model

    ### training ###
    num_epochs = 60
    batch_size = 16
    loss_fn = nn.BCELoss()    
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001, weight_decay=5e-4)
    train_x, val_x, train_y, val_y = train_test_split(input_data, truth_data, test_size=0.2, random_state=42)

    train_loader = DataLoader(list(zip(train_x, train_y)), batch_size=batch_size, shuffle=True, num_workers=4)
    val_loader = DataLoader(list(zip(val_x, val_y)), batch_size=batch_size, shuffle=False, num_workers=4)
    print('length of train_Set:', len(train_x))
    print('length of val_Set:', len(val_x))

    rmse_train_list, rmse_val_list, rp_train, rp_val = train_model(
        model, num_epochs, optimizer, train_loader, val_loader, loss_fn)

    ### save model ###
    torch.save(model, f'{result_dir}/model.pt')

    # PCC Plot
    plt.figure(figsize=(10, 6))
    plt.plot(rp_train, label='Training PCC')
    plt.plot(rp_val, label = 'Val PCC')
    plt.title('Pearson Correlation Coefficient over epochs')
    plt.xlabel('Epochs')
    plt.ylabel('PCC')
    plt.legend()
    plt.savefig(f'{result_dir}/pcc_plot.png')  # Save the plot as a PNG file

    # RMSE Plot
    plt.figure(figsize=(10, 6))
    plt.plot(rmse_train_list, label='train RMSE')
    plt.plot(rmse_val_list, label='val RMSE')
    plt.title('Root Mean Square Error over epochs')
    plt.xlabel('Epochs')
    plt.ylabel('RMSE')
    plt.legend()
    plt.savefig(f'{result_dir}/rmse_plot.png')  # Save the plot as a PNG file

    y_true, y_pred = test_basemodel(model, val_loader, device)
    plt.figure()
    plt.scatter(y_true, y_pred, alpha=0.5)
    plt.xlabel('Experimental HDX')
    plt.ylabel('Predicted HDX')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.plot([0, 1], [0, 1], '--', color='gray')
    plt.savefig(f'{result_dir}/val_result.png')