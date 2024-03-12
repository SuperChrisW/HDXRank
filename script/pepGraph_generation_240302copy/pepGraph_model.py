import torch
import torch.nn as nn
import torch.nn.functional as F
import torch_geometric as tg
from torch_geometric.nn import GCNConv, SAGEConv


class BiLSTM(nn.Module):
    def __init__(self, args):
        super(BiLSTM, self).__init__()
        # linear layer initialization
        self.seq_dim = args['seq_dim']
        self.struc_dim = args['struc_dim']
        self.evo_dim = args['evo_dim']
        self.time_dim = args['time_dim']

        #Convolutional layers
        self.hidden_channels = args['num_hidden_channels']
        self.out_channels = args['num_out_channels']
        self.drop_rate = args['drop_out']
        
        self.seq_linear = nn.Linear(10, self.seq_dim) # seq based feature: 10
        self.struc_linear = nn.Linear(5, self.struc_dim) # struc based feature: 5
        self.evo_linear = nn.Linear(30, self.evo_dim) # evo based feature: 30
        self.time_linear = nn.Linear(1, self.time_dim) # time based feature: 1

        self.conv1 = nn.Conv2d(in_channels=1, out_channels=self.hidden_channels, kernel_size=(5, 5), padding='same')
        self.bn1 = nn.BatchNorm2d(self.hidden_channels)
        self.conv2 = nn.Conv2d(self.hidden_channels, self.out_channels, kernel_size=(7, 7), padding='same')
        self.bn2 = nn.BatchNorm2d(self.out_channels)

        stride = 1
        kernel_size = 3
        input_width = self.seq_dim + self.struc_dim + self.evo_dim + self.time_dim
        input_height = 30
        output_height = (input_height - kernel_size) // stride + 1
        output_width = (input_width - kernel_size) // stride + 1
        self.pool = nn.MaxPool2d(kernel_size=(3, 3), stride=(1, 1), padding=0)
        self.dropout = nn.Dropout(self.drop_rate)

        self.lstm = nn.LSTM(input_size=output_height*self.out_channels, hidden_size=8, batch_first=True, bidirectional=True)
        self.fc1 = nn.Linear(16, 16)  # 16 units total (8 in each direction)

    def reset_parameters(self):
        self.seq_linear.reset_parameters()
        self.struc_linear.reset_parameters()
        self.evo_linear.reset_parameters()
        self.time_linear.reset_parameters()
        self.conv1.reset_parameters()
        self.conv2.reset_parameters()
        self.lstm.reset_parameters()
        self.fc1.reset_parameters()

    def forward(self, x):
        ### Featrue transformation from different sources: sequence, structure, evolution, time
        x_seq = F.relu(self.seq_linear(x[:, :, :, 35:45]))
        x_struc = F.relu(self.struc_linear(x[:, :, :, 0:5]))
        x_evo = F.relu(self.evo_linear(x[:, :, :, 5:35]))
        x_time = x[:, :, :, -1].unsqueeze(-1)
        x_time = F.relu(self.time_linear(x_time))

        x = torch.cat((x_seq, x_struc, x_evo, x_time), dim=3)
        x = x.permute(0, 1, 3, 2)
        ### Convolutional layers
        x = F.relu(self.bn1(self.conv1(x)))
        x = F.relu(self.bn2(self.conv2(x)))
        x = self.pool(x)

        x = x.permute(0, 2, 1, 3)  # Flatten for LSTM
        x = x.reshape(-1, x.shape[1], x.shape[2]*x.shape[3])

        ### LSTM layers
        x = self.dropout(x)
        x, _ = self.lstm(x)  # Reshape for LSTM
        x = self.dropout(x[:, -1, :])  # Use output of last LSTM sequence
        x = F.relu(self.fc1(x))
        return x

class GCN(nn.Module):
    def __init__(self, args):
        super().__init__()
        self.in_dim = args['in_dim']
        self.hidden_dim = args['GNN_hidden_dim']
        self.module_out_dim = args['GNN_out_dim']
        self.seq_dim = args['seq_dim']
        self.struc_dim = args['struc_dim']
        self.evo_dim = args['evo_dim']
        self.time_dim = args['time_dim']
        self.drop_rate = args['drop_out']
        self.num_GNN_layers = args['num_GNN_layers']

        gcn_layers = []
        GNN_layer = {'GCN': GCNConv, 'GraphSAGE': SAGEConv}
        if args['GNN_type'] == 'GCN':
            gcn_layers.append(GCNConv(self.in_dim-1, self.hidden_dim))
            for _ in range(self.num_GNN_layers - 2):
                gcn_layers.append(GCNConv(self.hidden_dim, self.hidden_dim))
            gcn_layers.append(GCNConv(self.hidden_dim, self.module_out_dim))
        elif args['GNN_type'] == 'GraphSAGE':
            gcn_layers.append(SAGEConv(self.in_dim-1, self.hidden_dim, aggr = 'mean'))
            for _ in range(self.num_GNN_layers - 2):
                gcn_layers.append(SAGEConv(self.hidden_dim, self.hidden_dim, aggr = 'mean'))
            gcn_layers.append(SAGEConv(self.hidden_dim, self.module_out_dim, aggr = 'mean'))
        self.gcn_layers = nn.ModuleList(gcn_layers)

        self.drop_out = nn.Dropout(self.drop_rate)
        self.bn = tg.nn.global_mean_pool
        self.fc = nn.Linear(self.module_out_dim, 1)

    def forward(self, graph_data):
        #graph_data.edge_attr = torch.tanh(-graph_data.edge_attr/2+2)+1 # transform distance to edge feature
        x = graph_data.x
        batch = graph_data.batch
        edge_index = graph_data.edge_index
        for layer in self.gcn_layers:
            x = F.relu(layer(x, edge_index))
            x = self.drop_out(x)
        x = self.bn(x, batch)
        x = F.sigmoid(self.fc(x))
        return x
    

class MixGCNBiLSTM(nn.Module):
    def __init__(self, args):
        super(MixGCNBiLSTM, self).__init__()
        self.in_dim = args['in_dim']
        self.hidden_dim = args['final_hidden_dim']
        self.module_out_dim = args['GNN_out_dim']
        self.drop_out = args['drop_out']

        self.bilstm = BiLSTM(args)
        self.gcn = GCN(args)

        self.fc1 = nn.Linear(self.module_out_dim*2, self.hidden_dim)
        self.fc2 = nn.Linear(self.hidden_dim, 1)
        self.dropout = nn.Dropout(self.drop_out)

    def forward(self, graph_data, seq_embedding):
        x_bilstm = self.bilstm(seq_embedding)
        x_gcn = self.gcn(graph_data)

        x_integrated = torch.cat((x_bilstm, x_gcn), dim=1)
        x = self.dropout(F.relu(self.fc1(x_integrated)))
        x = torch.sigmoid(self.fc2(x))
        return x