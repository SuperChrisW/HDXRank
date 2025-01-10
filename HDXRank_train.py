"""
2025/1/8
Author: WANG Liyao
Paper: HDXRank: A Deep Learning Framework for Ranking Protein complex predictions with Hydrogen Deuterium Exchange Data
Note: 
Modular training pipeline for GearNet model
"""
import os
import torch
import argparse
import logging
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from tqdm import tqdm
from torch import nn
from torchdrug.data import DataLoader
from GearNet import GearNet
from HDXRank_utilis import XML_process

# Utility functions
def load_data(tasks):
    """
    Load and preprocess data.

    Args:
        tasks (dict): Parsed XML file.

    Returns:
        tuple: apo_input, complex_input
    """
    summary_HDX_file = os.path.join(tasks["GeneralParameters"]["RootDir"], f"{tasks['GeneralParameters']['TaskFile']}.xlsx")
    hdx_df = pd.read_excel(summary_HDX_file, sheet_name='Sheet1')
    hdx_df = hdx_df.dropna(subset=['structure_file']).drop_duplicates(subset=['structure_file'])

    pepGraph_dir = tasks['GeneralParameters']['pepGraphDir']

    apo_input, complex_input = [], []

    logging.info('Loading data...')
    for _, row in tqdm(hdx_df.iterrows(), total=len(hdx_df)):
        pdb = row['structure_file'].strip().split('.')[0].upper()
        pepGraph_file = os.path.join(pepGraph_dir, f'{pdb}.pt')

        if os.path.isfile(pepGraph_file):
            pepGraph_ensemble = torch.load(pepGraph_file)
            if row['complex_state'] == 'single':
                apo_input.extend(pepGraph_ensemble)
            else:
                complex_input.extend(pepGraph_ensemble)

    logging.info(f"Length of apo data: {len(apo_input)}")
    logging.info(f"Length of complex data: {len(complex_input)}")

    return apo_input, complex_input

def prepare_model(input_dim, hidden_dims, num_relation, device):
    """
    Prepare the model and optimizer.

    Args:
        input_dim (int): Input dimension.
        hidden_dims (list): Hidden dimensions for the model.
        num_relation (int): Number of relations.
        device (torch.device): Device to use.

    Returns:
        tuple: model, optimizer, loss_fn
    """
    model = GearNet(input_dim=input_dim, hidden_dims=hidden_dims, num_relation=num_relation,
                    batch_norm=True, concat_hidden=True, readout='sum', activation='relu', short_cut=True).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001, weight_decay=5e-4)
    loss_fn = nn.BCELoss()
    return model, optimizer, loss_fn

def train_model(model, optimizer, loss_fn, train_loader, device, num_epochs):
    """
    Train the model.

    Args:
        model: PyTorch model instance.
        optimizer: Optimizer instance.
        loss_fn: Loss function instance.
        train_loader: DataLoader for training data.
        device: Device for computation.
        num_epochs: Number of epochs.

    Returns:
        tuple: rmse_train_list, rp_train
    """
    rp_train, rmse_train_list = [], []

    for epoch in range(num_epochs):
        model.train()
        list1_train, list2_train = [], []
        epoch_train_losses = []

        for graph_batch in train_loader:
            graph_batch = graph_batch.to(device)
            targets = graph_batch.y
            node_feat = graph_batch.residue_feature.float()
            outputs = model(graph_batch, node_feat)

            train_loss = loss_fn(outputs, targets)
            optimizer.zero_grad()
            train_loss.backward()
            optimizer.step()

            epoch_train_losses.append(train_loss.item())
            targets = targets.detach().cpu().numpy()
            outputs = outputs.detach().cpu().numpy()
            list1_train.extend(targets)
            list2_train.extend(outputs)

        epoch_train_loss = np.mean(epoch_train_losses)
        epoch_rp_train = np.corrcoef(list2_train, list1_train)[0, 1]
        rp_train.append(epoch_rp_train)

        y = np.array(list1_train).reshape(-1, 1)
        x = np.array(list2_train).reshape(-1, 1)
        epoch_train_rmse = np.sqrt(((y - x) ** 2).mean())
        rmse_train_list.append(epoch_train_rmse)

        logging.info(f'Epoch {epoch}: Loss {epoch_train_loss:.3f}, rho {epoch_rp_train:.3f}, RMSE {epoch_train_rmse:.3f}')

    return rmse_train_list, rp_train

def save_checkpoint(model, optimizer, epoch, file_path):
    checkpoint = {
        'epoch': epoch,
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict()
    }
    torch.save(checkpoint, file_path)

def main():
    parser = argparse.ArgumentParser(description='Train new HDXRank model.')
    parser.add_argument('-input', type=str, required=True, help='path to XML task file (require general parameters)')
    parser.add_argument('-save', type=str, required=True, help='path to save the model')
    parser.add_argument('-epoch', type=int, default=100, help='Number of epochs to train')
    parser.add_argument('-cuda', type=int, default=0, help='CUDA device number')
    parser.add_argument('-train_val_split', type=float, default=0.2, help='Proportion of data to use for validation')
    parser.add_argument('-random_state', type=int, default=42, help='Random state for data splitting')
    parser.add_argument('-batch_size', type=int, default=16, help='Batch size for training')
    parser.add_argument('-repeat', type=int, default=1, help='repeat training process')
    args = parser.parse_args()

    device = torch.device(f'cuda:{args.cuda}' if torch.cuda.is_available() else 'cpu')

    tasks = XML_process(args.input)
    apo_input, complex_input = load_data(tasks)

    for i in range(args.repeat):
        train_apo, val_apo = train_test_split(apo_input, test_size=args.train_val_split, random_state=args.random_state)
        train_complex, val_complex = train_test_split(complex_input, test_size=args.train_val_split, random_state=args.random_state)

        train_set = train_apo + train_complex
        val_set = val_apo + val_complex
        train_loader = DataLoader(train_set+val_set, batch_size=args.batch_size, shuffle=True) # use all data to train the final model
        #val_loader = DataLoader(val_set, batch_size=args.batch_size, shuffle=False)
        print('training set:', len(train_set+val_set))
        #print('val set:', len(val_set))

        model, optimizer, loss_fn = prepare_model(input_dim=56, hidden_dims=[512, 512, 512], num_relation=7, device=device)
        rmse_train_list, rp_train = train_model(model, optimizer, loss_fn, train_loader, device, args.epoch)

        model_name = f'HDXRank_GN56_epoch{args.epoch}_v{i}.pth'
        save_checkpoint(model, optimizer, args.epoch, os.path.join(args.save, model_name))

if __name__ == "__main__":
    # Logging setup
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler()]
    )
    main()
