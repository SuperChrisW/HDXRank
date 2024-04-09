import torch
import os
import accelerate
from accelerate import Accelerator
from pepGraph_BiLSTM import GNN_v2

# load torch model using accelerate unwrap

model_path = '/home/lwang/models/HDX_LSTM/results/240402_GAT_v2/hop1/best_model.pt'
accelerator = Accelerator()

config = {
            'num_epochs': 300,
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
        'result_dir': '/home/lwang/models/HDX_LSTM/results/240402_GAT_v2/hop1'
}
model = GNN_v2(training_args)
unwrapped_model = accelerator.unwrap_model(model)
path_to_checkpoint = os.path.join(model_path,"pytorch_model.bin")
unwrapped_model.load_state_dict(torch.load(path_to_checkpoint))

#show the model architecture
print(unwrapped_model)