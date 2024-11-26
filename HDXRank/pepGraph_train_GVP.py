### training ###
import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

from comet_ml import Experiment
from comet_ml.integration.pytorch import log_model

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.loader import DataLoader
from GVP_model import MQAModel

import pandas as pd
import numpy as np
from tqdm import tqdm
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.metrics import r2_score, mean_squared_error
from scipy.stats import spearmanr

def train_model(model, num_epochs, optimizer, train_loader, val_loader, loss_fn, device, experiment, result_fpath, data_log = True):
    for epoch in range(num_epochs):
        list1_train = []
        list2_train = []
        epoch_train_losses = []

        model.train()
        for graph_batch in train_loader:
            optimizer.zero_grad()
            
            graph_batch = graph_batch.to(device)
            targets = graph_batch['residue'].y.to(torch.float32)

            nodes = (graph_batch['residue'].node_s, graph_batch['residue'].node_v)
            edges = (graph_batch['residue'].edge_s, graph_batch['residue'].edge_v)
            outputs = model(nodes, graph_batch.edge_index_dict, edges,
                            batch = graph_batch['residue'].batch)

            train_loss = loss_fn(outputs, targets)
            train_loss.backward()
            optimizer.step()

            epoch_train_losses.append(train_loss.item())
            targets = targets.detach().cpu()
            outputs = outputs.detach().cpu()
            list1_train=np.append(list1_train,targets)
            list2_train=np.append(list2_train,outputs)

        epoch_train_losses = np.mean(epoch_train_losses)
        x = np.array(list1_train).reshape(-1,1)
        y = np.array(list2_train).reshape(-1,1)
        epoch_rp_train = spearmanr(x, y)[0]
        epoch_r2_train = r2_score(x, y)       
        epoch_train_rmse= np.sqrt(((y - x) ** 2).mean())
        
        if data_log:
            experiment.log_metric('train_loss', epoch_train_losses, step = epoch)
            experiment.log_metric('train_scc', epoch_rp_train, step = epoch)
            experiment.log_metric('train_rmse', epoch_train_rmse, step = epoch)

        print('Epoch  Train_Loss  Train_SPR    Train_RMSE   Train_R2')
        print("{:5d}  {:10.3f}  {:9.3f}  {:10.3f}  {:10.3f}".format(
        epoch, epoch_train_losses, epoch_rp_train, epoch_train_rmse, epoch_r2_train))

        ### validation and early-stopping
        '''if epoch % 5 == 0:
            model.eval()
            epoch_val_losses = []
            val_true = []
            val_pred = []
            with torch.no_grad():
                for graph_batch in val_loader:
                    graph_batch = graph_batch.to(device)
                    targets = graph_batch['residue'].y.to(torch.float32)

                    nodes = (graph_batch['residue'].node_s, graph_batch['residue'].node_v)
                    edges = (graph_batch['residue'].edge_s, graph_batch['residue'].edge_v)
                    outputs = model(nodes, graph_batch.edge_index_dict, edges,
                                    batch = graph_batch['residue'].batch)
                
                    val_loss = loss_fn(outputs, targets)
                    val_true.extend(targets.detach().cpu().numpy())
                    val_pred.extend(outputs.detach().cpu().numpy())
                    epoch_val_losses.append(val_loss.item())
                val_losses_mean = np.mean(epoch_val_losses)

            if data_log:
                val_scc = spearmanr(val_true, val_pred)[0]
                val_rmse = np.sqrt(mean_squared_error(val_true, val_pred))
                val_r2 = r2_score(val_true, val_pred)
                experiment.log_metric('val_loss', val_losses_mean, step = epoch)
                experiment.log_metric('val_scc', val_scc, step = epoch)
                experiment.log_metric('val_rmse', val_rmse, step = epoch)
                experiment.log_metric('val_r2', val_r2, step = epoch)

        print('Epoch  Train_Loss  Train_SCC    Train_RMSE   Val_SCC    Val_RMSE    Val_R2')
        print("{:5d}  {:10.3f}  {:9.3f}  {:10.3f}  {:10.3f}  {:10.3f}  {:10.3f}".format(
        epoch, epoch_train_losses, epoch_rp_train, epoch_train_rmse, val_scc, val_rmse, val_r2))'''
        #if epoch % 10 == 0:
        #    save_checkpoint(model, optimizer, epoch, f'{result_fpath}_epoch{epoch}.pth')

    model.eval()
    complex_labels = []
    val_preds = []
    val_targets = []
    with torch.no_grad():
        for graph_batch in val_loader: 
            graph_batch = graph_batch.to(device)
            targets = graph_batch['residue'].y.to(torch.float32)

            nodes = (graph_batch['residue'].node_s, graph_batch['residue'].node_v)
            edges = (graph_batch['residue'].edge_s, graph_batch['residue'].edge_v)
            outputs = model(nodes, graph_batch.edge_index_dict, edges,
                            batch = graph_batch['residue'].batch)
            
            val_preds.extend(outputs.cpu().numpy())
            val_targets.extend(targets.cpu().numpy())
            complex_labels.extend(graph_batch['residue'].is_complex.cpu().numpy())
        
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
    single_r2 = r2_score(val_targets[single_mask], val_preds[single_mask])
    complex_r2 = r2_score(val_targets[complex_mask], val_preds[complex_mask])
    total_r2 = r2_score(val_targets, val_preds)
    print('------------ Model Validation ------------')
    print('Single RMSE   Single SPR   Single R2')
    print("{:10.3f}  {:9.3f}  {:10.3f}".format(single_rmse, single_spr, single_r2))
    print('Complex RMSE  Complex SPR  Complex R2')
    print("{:10.3f}  {:9.3f}  {:10.3f}".format(complex_rmse, complex_spr, complex_r2))
    print('Total RMSE    Total SPR    Total R2')
    print("{:10.3f}  {:9.3f}  {:10.3f}".format(total_rmse, total_spr, total_r2))
    
    save_checkpoint(model, optimizer, epoch, f'{result_fpath}_epoch{epoch}.pth')

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

    pepGraph_dir = os.path.join(root_dir, 'graph_ensemble_simpleGVP', training_args['cluster'])
    result_dir = training_args['result_dir']
    if not os.path.exists(result_dir):
        os.mkdir(result_dir)

    device = training_args['device']
    ### data preparation ###
    apo_input = []
    complex_input = []
    hdx_df = hdx_df.drop_duplicates(subset=['structure_file'])

    print('data loading')
    for row_id, row in tqdm(hdx_df.iterrows()):
        pdb = row['structure_file'].strip().split('.')[0].upper()
        pepGraph_file = f'{pepGraph_dir}/{pdb}.pt'
        if os.path.isfile(pepGraph_file):
            pepGraph_ensemble = torch.load(pepGraph_file)
            if row['complex_state'] == 'single':
                apo_input.extend(pepGraph_ensemble)
            else:
                complex_input.extend(pepGraph_ensemble)
        else:
            continue

    #apo_input = torch.load(f'{pepGraph_dir}/train_val_apo.pt')
    #complex_input = torch.load(f'{pepGraph_dir}/train_val_complex.pt')
    print('length of apo data:', len(apo_input))
    print('length of complex data:', len(complex_input))

    ### model initialization ###
    torch.cuda.empty_cache()

    #GVP-GNN outputing a scaler value for each graph
    node_in_dim = (56,3) # 78 or 98 or 1280
    node_h_dim = (1280, 12) # 392 or 1280
    edge_in_dim = (32,1)
    edge_h_dim = (128, 4)
    #metadata = apo_input[0][0].metadata()
    metadata = apo_input[0].metadata()
    remove_list = ['knn_edge', 'forward_1_edge', 'forward_2_edge', 'backward_1_edge', 'backward_2_edge']
    for edge in remove_list:
        metadata[1].remove(('residue', edge, 'residue')) if ('residue', edge, 'residue') in metadata[1] else None

    ### training ###
    '''dataset = apo_input + complex_input
    type_label = [0]*len(apo_input) + [1]*len(complex_input)
    loss_fn = nn.BCELoss()
    kf = StratifiedKFold(n_splits=config['cross_validation_num'], shuffle=True, random_state=42)

    for fold, (train_index, val_index) in enumerate(kf.split(dataset, type_label)):
        train_set = [dataset[i] for i in train_index]
        val_set = [dataset[i] for i in val_index]'''

    loss_fn = nn.BCELoss()
    for fold in range(config['cross_validation_num']):
        model = MQAModel(node_in_dim, node_h_dim, 
                 edge_in_dim, edge_h_dim, metadata,
                num_layers=config['num_GNN_layers'], drop_rate=training_args['drop_out']).to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=config['learning_rate'], weight_decay=config['weight_decay'])

        print(f'-----------------------------------------')
        print(f'-              Fold {fold+1}               -')
        print(f'-----------------------------------------')

        train_apo, val_apo = train_test_split(apo_input, test_size=0.2, random_state=42)
        train_complex, val_complex = train_test_split(complex_input, test_size=0.2, random_state=42)

        train_set = train_apo + train_complex
        val_set = val_apo + val_complex

        train_loader = DataLoader(train_set, batch_size = config['batch_size'], shuffle=True, num_workers=config['num_workers'])
        val_loader =  DataLoader(val_set, batch_size = config['batch_size'], shuffle=False, num_workers=config['num_workers'])
        print('length of train_Set:', len(train_set))
        print('length of val_Set:', len(val_set))

        # train and save model at checkpoints
        fname = training_args['file_name']+'_v'+str(fold)
        result_fpath = os.path.join(training_args['result_dir'], fname)
        train_model(
            model, config['num_epochs'], optimizer, train_loader, val_loader, loss_fn, device,
            experiment, result_fpath, training_args['data_log'])
        if training_args['data_log']:
            log_model(experiment, model=model, model_name = 'PEP-HDX')
        
        break

if __name__ == "__main__":
    cluster = 'cluster1_8A_manual_rescale'
    device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')
    config = {
            'num_epochs':100,
            'batch_size': 16,
            'learning_rate': 0.001,
            'weight_decay': 5e-4,
            'GNN_type': f'model_GVP56_{cluster}_layer3-knnE',
            'num_GNN_layers': 3,
            'cross_validation_num': 1,
            'num_workers': 4,
    }

    training_args = {'num_hidden_channels': 10, 'num_out_channels': 20, 
            'feat_in_dim': 56, 'topo_in_dim': 0, 'num_heads': 8, 'GNN_hidden_dim': 32,
            'GNN_out_dim': 16, 'LSTM_out_dim': 16, 'num_GNN_layers': config['num_GNN_layers'],

            'final_hidden_dim': 16,

            'drop_out': 0.5, 'num_GNN_layers': config['num_GNN_layers'], 'GNN_type': config['GNN_type'],
            'graph_hop': 'hop1', 'batch_size': config['batch_size'],
            'result_dir': '/home/lwang/models/HDX_LSTM/results/240817_GVP',
            'file_name': f'model_GVP56_{cluster}_layer3-knnE',
            'data_log': False,
            'device': device,
            'cluster': cluster
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