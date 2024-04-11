import torch
import torch.nn as nn
import torch.nn.functional as F
import torch_geometric as tg
from torch_geometric.nn import GCNConv, SAGEConv, GATConv, global_mean_pool as gmp


class BiLSTM(nn.Module):
    def __init__(self, args):
        super(BiLSTM, self).__init__()
        #Convolutional layers
        self.hidden_channels = args['num_hidden_channels']
        self.out_channels = args['num_out_channels']
        self.module_out_dim = args['LSTM_out_dim']
        self.drop_rate = args['drop_out']

        self.conv1 = nn.Conv2d(in_channels=1, out_channels=self.hidden_channels, kernel_size=(5, 5), padding='same')
        self.bn1 = nn.BatchNorm2d(self.hidden_channels)
        self.conv2 = nn.Conv2d(self.hidden_channels, self.out_channels, kernel_size=(7, 7), padding='same')
        self.bn2 = nn.BatchNorm2d(self.out_channels)

        stride = 1
        kernel_size = 3
        input_width = args['feat_in_dim'] +2 # two extra feat: log_t and sequence length/ max_len
        input_height = 30
        output_height = (input_height - kernel_size) // stride + 1
        output_width = (input_width - kernel_size) // stride + 1
        self.pool = nn.MaxPool2d(kernel_size=(3, 3), stride=(1, 1), padding=0)
        self.dropout = nn.Dropout(self.drop_rate)

        self.lstm = nn.LSTM(input_size=output_height*self.out_channels, hidden_size=self.module_out_dim, batch_first=True, bidirectional=True)
        self.fc1 = nn.Linear(self.module_out_dim * 2, self.module_out_dim)  # 16 units total (8 in each direction)

    def reset_parameters(self):
        self.conv1.reset_parameters()
        self.conv2.reset_parameters()
        self.lstm.reset_parameters()
        self.fc1.reset_parameters()

    def forward(self, x):
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

class GNN_v2(nn.Module):
    def __init__(self, args):
        super().__init__()
        self.feat_in_dim = args['feat_in_dim']
        self.topo_in_dim = args['topo_in_dim']
        self.hidden_dim = args['GNN_hidden_dim']
        self.module_out_dim = args['GNN_out_dim']
        self.num_heads = args['num_heads']
        self.drop_rate = args['drop_out']
        self.num_GNN_layers = args['num_GNN_layers']

        # Separate GAT layers for feature and topology matrices
        self.gat_layers_feat = nn.ModuleList([
            GATConv(self.feat_in_dim, self.hidden_dim, heads=self.num_heads, dropout=self.drop_rate) if i == 0 else
            GATConv(self.hidden_dim * self.num_heads, self.hidden_dim, heads=self.num_heads) if i < self.num_GNN_layers - 1 else
            GATConv(self.hidden_dim * self.num_heads, self.module_out_dim, heads=self.num_heads)
            for i in range(self.num_GNN_layers)
        ])

        self.gat_layers_topo = nn.ModuleList([
            GATConv(self.topo_in_dim, self.hidden_dim, heads=self.num_heads, dropout=self.drop_rate) if i == 0 else
            GATConv(self.hidden_dim * self.num_heads, self.hidden_dim, heads=self.num_heads) if i < self.num_GNN_layers - 1 else
            GATConv(self.hidden_dim * self.num_heads, self.module_out_dim, heads=self.num_heads)
            for i in range(self.num_GNN_layers)
        ])

        self.drop_out = nn.Dropout(self.drop_rate)
        self.fc = nn.Linear(self.module_out_dim * self.num_heads * 2, self.module_out_dim)  # Adjusted for concatenated vector

    def reset_parameters(self):
        for layer in self.gat_layers_feat:
            layer.reset_parameters()
        for layer in self.gat_layers_topo:
            layer.reset_parameters()
        self.fc.reset_parameters()

    def forward_branch(self, branch_layers, x, edge_index):
        for layer in branch_layers:
            x = F.relu(layer(x, edge_index))
            x = self.drop_out(x)
        return x

    def forward(self, graph_data):
        topo_x = graph_data.x[:, self.feat_in_dim:self.feat_in_dim + self.topo_in_dim]
        feat_x = graph_data.x[:, :self.feat_in_dim]
        
        batch = graph_data.batch
        edge_index = graph_data.edge_index

        # Process feature and topology matrices separately
        processed_feat = self.forward_branch(self.gat_layers_feat, feat_x, edge_index)
        processed_topo = self.forward_branch(self.gat_layers_topo, topo_x, edge_index)

        # Concatenate processed features and topology matrices
        x = torch.cat((processed_feat, processed_topo), dim=-1)

        # Global mean pooling
        x = gmp(x, batch)
        
        # Fully connected layer
        x = torch.sigmoid(self.fc(x))
        
        return x
    
class MixGCNBiLSTM(nn.Module):
    def __init__(self, args):
        super(MixGCNBiLSTM, self).__init__()
        self.feat_in_dim = args['feat_in_dim']
        self.topo_in_dim = args['topo_in_dim']
        self.hidden_dim = args['final_hidden_dim']
        self.module_out_dim = args['GNN_out_dim'] + args['LSTM_out_dim']
        self.drop_out = args['drop_out']
        self.batch_size = args['batch_size']

        self.bilstm = BiLSTM(args)
        self.gnn = GNN_v2(args)

        self.fc1 = nn.Linear(self.module_out_dim, self.hidden_dim)
        self.fc2 = nn.Linear(self.hidden_dim, 1)
        self.dropout = nn.Dropout(self.drop_out)
    
    def reset_parameters(self):
        self.bilstm.reset_parameters()
        self.gnn.reset_parameters()
        self.fc1.reset_parameters()
        self.fc2.reset_parameters()

    def forward(self, graph_data):
        seq_embedding = graph_data['seq_embedding'].reshape(-1, 1, 30, self.feat_in_dim+2)
        x_bilstm = self.bilstm(seq_embedding)
        x_gcn = self.gnn(graph_data)
        x_integrated = torch.cat((x_bilstm, x_gcn), dim=1)
        x = self.dropout(F.relu(self.fc1(x_integrated)))
        x = torch.sigmoid(self.fc2(x))
        return x.squeeze(-1)