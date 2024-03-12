import torch
import numpy as np
from pepGraph_model import MixGCNBiLSTM, BiLSTM, GCN

args = {'in_dim': 46, 'final_hidden_dim': 48,
    'seq_dim': 10, 'struc_dim': 10, 'evo_dim': 10, 'time_dim': 5,
    'num_hidden_channels': 10, 'num_out_channels': 20, 
    'GNN_out_dim': 16, 'GNN_hidden_dim': 16,
    'drop_out': 0.5, 'num_GNN_layers': 3, 'GNN_type': 'GCN'}

model = BiLSTM(args)
input  = torch.randn(2, 1, 30, 46)
print(input.shape)
output = model(input)