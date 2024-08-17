import torch
import torch.nn as nn
import torch.nn.functional as F
import torch_geometric as tg
from torch_geometric.nn import GCNConv, SAGEConv, GATConv, global_add_pool as gmp


class BiLSTM(nn.Module):
    def __init__(self, args):
        super(BiLSTM, self).__init__()

        #Convolutional layers
        self.hidden_channels = args['num_hidden_channels']
        self.out_channels = args['num_out_channels']
        self.drop_rate = args['drop_out']

        self.conv1 = nn.Conv2d(in_channels=1, out_channels=self.hidden_channels, kernel_size=(3, 3), padding='same')
        self.bn1 = nn.BatchNorm2d(self.hidden_channels)
        self.conv2 = nn.Conv2d(self.hidden_channels, self.out_channels, kernel_size=(5, 5), padding='same')
        self.bn2 = nn.BatchNorm2d(self.out_channels)

        stride = 1
        kernel_size = 3
        input_width = 48
        input_height = 30
        output_height = (input_height - kernel_size) // stride + 1
        output_width = (input_width - kernel_size) // stride + 1
        self.pool = nn.MaxPool2d(kernel_size=(3, 3), stride=(1, 1), padding=0)
        self.dropout = nn.Dropout(self.drop_rate)

        self.lstm = nn.LSTM(input_size=output_height*self.out_channels, hidden_size=8, batch_first=True, bidirectional=True)
        self.fc1 = nn.Linear(16, 10)  # 16 units total (8 in each direction)
        self.fc2 = nn.Linear(10, 1)  # 16 units total (8 in each direction)

    def reset_parameters(self):
        self.conv1.reset_parameters()
        self.conv2.reset_parameters()
        self.lstm.reset_parameters()
        self.fc1.reset_parameters()

    def forward(self, graph):
        ### Featrue transformation from different sources: sequence, structure, evolution, time
        x = graph['seq_embedding'].reshape(-1, 1, 30, 48)
        '''
        SASA = x[:, :, :, 0:1]
        HHblites = x[:, :, :, 5:35]
        HDMD = x[:, :, :, 35:40]
        x = torch.cat((SASA, HDMD, HHblites), dim=3)
        x = x.reshape(-1, 1, 30, 36)
        '''
        #x = x.permute(0, 1, 3, 2)
        ### Convolutional layers
        x = F.relu(self.bn1(self.conv1(x)))
        x = F.relu(self.bn2(self.conv2(x)))
        x = self.pool(x)
        x = x.permute(0, 3, 1, 2)  # Flatten for LSTM
        x = x.reshape(-1, x.shape[1], x.shape[2]*x.shape[3])
        ### LSTM layers
        x = self.dropout(x)
        x, _ = self.lstm(x) # Reshape for LSTM
        x = self.dropout(x[:, -1, :])  # Use output of last LSTM sequence
        x = F.relu(self.fc1(x))
        x = F.sigmoid(self.fc2(x))
        return x.squeeze(-1)

class BiLSTM_36(nn.Module):
    def __init__(self, args):
        super(BiLSTM, self).__init__()

        #Convolutional layers
        self.hidden_channels = args['num_hidden_channels']
        self.out_channels = args['num_out_channels']
        self.drop_rate = args['drop_out']

        self.conv1 = nn.Conv2d(in_channels=1, out_channels=self.hidden_channels, kernel_size=(3, 3), padding='same')
        self.bn1 = nn.BatchNorm2d(self.hidden_channels)
        self.conv2 = nn.Conv2d(self.hidden_channels, self.out_channels, kernel_size=(5, 5), padding='same')
        self.bn2 = nn.BatchNorm2d(self.out_channels)

        stride = 1
        kernel_size = 3
        input_width = 36
        input_height = 30
        output_height = (input_height - kernel_size) // stride + 1
        output_width = (input_width - kernel_size) // stride + 1
        self.pool = nn.MaxPool2d(kernel_size=(3, 3), stride=(1, 1), padding=0)
        self.dropout = nn.Dropout(self.drop_rate)

        self.lstm = nn.LSTM(input_size=output_height*self.out_channels, hidden_size=8, batch_first=True, bidirectional=True)
        self.fc1 = nn.Linear(16, 10)  # 16 units total (8 in each direction)
        self.fc2 = nn.Linear(10, 1)  # 16 units total (8 in each direction)

    def reset_parameters(self):
        self.conv1.reset_parameters()
        self.conv2.reset_parameters()
        self.lstm.reset_parameters()
        self.fc1.reset_parameters()

    def forward(self, graph):
        ### Featrue transformation from different sources: sequence, structure, evolution, time
        x = graph['seq_embedding'].reshape(-1, 1, 30, 48)

        SASA = x[:, :, :, 0:1]
        HHblites = x[:, :, :, 5:35]
        HDMD = x[:, :, :, 35:40]
        x = torch.cat((SASA, HDMD, HHblites), dim=3)
        x = x.reshape(-1, 1, 30, 36)

        x = x.permute(0, 1, 3, 2)
        ### Convolutional layers
        x = F.relu(self.bn1(self.conv1(x)))
        x = F.relu(self.bn2(self.conv2(x)))
        x = self.pool(x)
        x = x.permute(0, 2, 1, 3)  # Flatten for LSTM
        x = x.reshape(-1, x.shape[1], x.shape[2]*x.shape[3])
        ### LSTM layers
        x = self.dropout(x)
        x, _ = self.lstm(x) # Reshape for LSTM
        x = self.dropout(x[:, -1, :])  # Use output of last LSTM sequence
        x = F.relu(self.fc1(x))
        x = F.sigmoid(self.fc2(x))
        return x.squeeze(-1)

class GAT(nn.Module):
    def __init__(self, args):
        super().__init__()
        self.feat_in_dim = args['feat_in_dim']
        self.hidden_dim = args['GNN_hidden_dim']
        self.module_out_dim = args['GNN_out_dim']
        self.num_heads = args['num_heads']
        self.drop_rate = args['drop_out']
        self.num_GNN_layers = args['num_GNN_layers']

        # GAT layers
        self.gat_layers = nn.ModuleList([
            GATConv(self.feat_in_dim, self.hidden_dim, heads=self.num_heads, dropout=self.drop_rate) if i == 0 else
            GATConv(self.hidden_dim * self.num_heads, self.hidden_dim, heads=self.num_heads, dropout=self.drop_rate) if i < self.num_GNN_layers - 1 else
            GATConv(self.hidden_dim * self.num_heads, self.module_out_dim, heads=self.num_heads, dropout=self.drop_rate)
            for i in range(self.num_GNN_layers)
        ])

        self.batch_norms = nn.ModuleList([
            nn.BatchNorm1d(self.hidden_dim * self.num_heads) for _ in range(self.num_GNN_layers - 1)
        ])
        
        self.drop_out = nn.Dropout(self.drop_rate)
        self.fc = nn.Linear(self.module_out_dim * self.num_heads, 1)

    def forward(self, graph_data):
        x = graph_data['residue'].node_s['residue'][:, :self.feat_in_dim]
        num_nodes = graph_data['residue'].num_nodes
        self_loops = torch.arange(0,num_nodes, device = x.device)
        self_loops = self_loops.unsqueeze(0).repeat(2,1)
        edge_index = graph_data.edge_index_dict[('residue', 'radius_edge', 'residue')]
        edge_index = torch.cat([edge_index, self_loops], dim=1)
        batch = graph_data['residue'].batch

        for i, layer in enumerate(self.gat_layers):
            x = F.relu(layer(x, edge_index))
            if i < len(self.batch_norms):
                x = self.batch_norms[i](x)
            x = self.drop_out(x)

        x = gmp(x, batch)
        x = self.fc(x)
        x = torch.sigmoid(x)
        
        return x.squeeze(-1)
