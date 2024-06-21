import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_scatter import scatter_add
from torchdrug import layers


class BiLSTM(nn.Module):
    def __init__(self, args):
        super(BiLSTM, self).__init__()
        #Convolutional layers
        self.hidden_channels = args['num_hidden_channels']
        self.out_channels = args['num_out_channels']
        self.module_out_dim = args['LSTM_out_dim']
        self.drop_rate = args['drop_out']
        self.feat_in_dim = args['feat_in_dim']

        self.conv1 = nn.Conv2d(in_channels=1, out_channels=self.hidden_channels, kernel_size=(5, 5), padding='same')
        self.bn1 = nn.BatchNorm2d(self.hidden_channels)
        self.conv2 = nn.Conv2d(self.hidden_channels, self.out_channels, kernel_size=(7, 7), padding='same')
        self.bn2 = nn.BatchNorm2d(self.out_channels)

        stride = 1
        kernel_size = 3
        input_width = args['feat_in_dim'] +1 # one extra feat: sequence length/ max_len
        input_height = 30
        output_height = (input_height - kernel_size) // stride + 1
        output_width = (input_width - kernel_size) // stride + 1
        self.pool = nn.MaxPool2d(kernel_size=(3, 3), stride=(1, 1), padding=0)
        self.dropout = nn.Dropout(self.drop_rate)

        self.lstm = nn.LSTM(input_size=output_height*self.out_channels, hidden_size=self.module_out_dim, batch_first=True, bidirectional=True)
        self.fc1 = nn.Linear(self.module_out_dim * 2, self.module_out_dim)  # 16 units total (8 in each direction)
        self.output_fc = nn.Linear(self.module_out_dim, 1)

    def reset_parameters(self):
        self.conv1.reset_parameters()
        self.conv2.reset_parameters()
        self.lstm.reset_parameters()
        self.fc1.reset_parameters()
        self.output_fc.reset_parameters()

    def forward(self, x):
        #x = x.reshape(-1, 1, 30, self.feat_in_dim+1)
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
        x = F.sigmoid(self.output_fc(x)).squeeze(-1) # use if the model only has bilstm module
        return x

class GearNet(nn.Module):
    """
    Geometry Aware Relational Graph Neural Network proposed in
    `Protein Representation Learning by Geometric Structure Pretraining`_.

    .. _Protein Representation Learning by Geometric Structure Pretraining:
        https://arxiv.org/pdf/2203.06125.pdf

    Parameters:
        input_dim (int): input dimension
        hidden_dims (list of int): hidden dimensions
        num_relation (int): number of relations
        edge_input_dim (int, optional): dimension of edge features
        num_angle_bin (int, optional): number of bins to discretize angles between edges.
            The discretized angles are used as relations in edge message passing.
            If not provided, edge message passing is disabled.
        short_cut (bool, optional): use short cut or not
        batch_norm (bool, optional): apply batch normalization or not
        activation (str or function, optional): activation function
        concat_hidden (bool, optional): concat hidden representations from all layers as output
        readout (str, optional): readout function. Available functions are ``sum`` and ``mean``.
    """

    def __init__(self, input_dim, hidden_dims, num_relation, edge_input_dim=None, num_angle_bin=None,
                 short_cut=False, batch_norm=False, activation="relu", concat_hidden=False, readout="sum", args = None):
        super(GearNet, self).__init__()

        self.input_dim = input_dim
        self.output_dim = sum(hidden_dims) if concat_hidden else hidden_dims[-1]
        self.dims = [input_dim] + list(hidden_dims)
        self.edge_dims = [edge_input_dim] + self.dims[:-1]
        self.num_relation = num_relation
        self.num_angle_bin = num_angle_bin
        self.short_cut = short_cut
        self.concat_hidden = concat_hidden
        self.batch_norm = batch_norm
        self.module_out_dim = args['GNN_out_dim'] if args is not None else 64

        self.layers = nn.ModuleList()
        for i in range(len(self.dims) - 1):
            self.layers.append(layers.GeometricRelationalGraphConv(self.dims[i], self.dims[i + 1], num_relation,
                                                                   None, batch_norm, activation))
        if num_angle_bin:
            self.spatial_line_graph = layers.SpatialLineGraph(num_angle_bin)
            self.edge_layers = nn.ModuleList()
            for i in range(len(self.edge_dims) - 1):
                self.edge_layers.append(layers.GeometricRelationalGraphConv(
                    self.edge_dims[i], self.edge_dims[i + 1], num_angle_bin, None, batch_norm, activation))

        if batch_norm:
            self.batch_norms = nn.ModuleList()
            for i in range(len(self.dims) - 1):
                self.batch_norms.append(nn.BatchNorm1d(self.dims[i + 1]))

        if readout == "sum":
            self.readout = layers.SumReadout()
        elif readout == "mean":
            self.readout = layers.MeanReadout()
        else:
            raise ValueError("Unknown readout `%s`" % readout)

        # MLP output layer
        if concat_hidden:
            self.mlp = layers.MLP(sum(hidden_dims), self.module_out_dim,
                    batch_norm=False, dropout=0)
        else:
            self.mlp = layers.MLP(hidden_dims[-1], self.module_out_dim,
                        batch_norm=False, dropout=0, activation = 'relu')

    def forward(self, graph, input, all_loss=None, metric=None):
        """
        Compute the node representations and the graph representation(s).

        Parameters:
            graph (Graph): :math:`n` graph(s)
            input (Tensor): input node representations
            all_loss (Tensor, optional): if specified, add loss to this tensor
            metric (dict, optional): if specified, output metrics to this dict

        Returns:
            dict with ``node_feature`` and ``graph_feature`` fields:
                node representations of shape :math:`(|V|, d)`, graph representations of shape :math:`(n, d)`
        """
        hiddens = []
        layer_input = input
        if self.num_angle_bin:
            line_graph = self.spatial_line_graph(graph)
            edge_input = line_graph.node_feature.float()

        for i in range(len(self.layers)):
            hidden = self.layers[i](graph, layer_input)
            if self.short_cut and hidden.shape == layer_input.shape:
                hidden = hidden + layer_input
            if self.num_angle_bin:
                edge_hidden = self.edge_layers[i](line_graph, edge_input)
                edge_weight = graph.edge_weight.unsqueeze(-1)
                node_out = graph.edge_list[:, 1] * self.num_relation + graph.edge_list[:, 2]
                update = scatter_add(edge_hidden * edge_weight, node_out, dim=0,
                                     dim_size=graph.num_node * self.num_relation)
                update = update.view(graph.num_node, self.num_relation * edge_hidden.shape[1])
                update = self.layers[i].linear(update)
                update = self.layers[i].activation(update)
                hidden = hidden + update
                edge_input = edge_hidden
            if self.batch_norm:
                hidden = self.batch_norms[i](hidden)
            hiddens.append(hidden)
            layer_input = hidden

        if self.concat_hidden:
            node_feature = torch.cat(hiddens, dim=-1)
        else:
            node_feature = hiddens[-1]
        graph_feature = self.readout(graph, node_feature)
        pred = self.mlp(graph_feature).squeeze(-1)

        return F.relu(pred)
    
class MixBiLSTM_GearNet(nn.Module):
    def __init__(self, args):
        super(MixBiLSTM_GearNet, self).__init__()
        self.feat_in_dim = args['feat_in_dim']
        self.topo_in_dim = args['topo_in_dim']
        self.hidden_dim = args['final_hidden_dim']
        self.module_out_dim = args['GNN_out_dim'] + args['LSTM_out_dim']
        self.drop_out = args['drop_out']
        self.batch_size = args['batch_size']

        self.bilstm = BiLSTM(args)
        #self.gnn = GearNet(input_dim = self.feat_in_dim+self.topo_in_dim, hidden_dims = [512,512,512],
        #            num_relation=7, batch_norm=True, concat_hidden=True, readout='sum', activation = 'relu', short_cut=True, args = args)

        self.gnn = GearNet(input_dim=self.feat_in_dim+self.topo_in_dim, hidden_dims=[512, 512, 512], 
                              num_relation=7, edge_input_dim=59, num_angle_bin=8,
                              batch_norm=True, concat_hidden=True, short_cut=True, readout="sum", activation = 'relu')

        self.attention_fc = nn.Linear(1, 16)
        self.fc1 = nn.Linear(16, self.hidden_dim)
        self.fc2 = nn.Linear(self.hidden_dim, 1)
        self.dropout = nn.Dropout(self.drop_out)
    
    def reset_parameters(self):
        self.bilstm.reset_parameters()
        self.gnn.reset_parameters()
        self.fc1.reset_parameters()
        self.fc2.reset_parameters()

    def forward(self, graph_data):
        seq_embedding = graph_data.seq_embedding.reshape(-1, 1, 30, self.feat_in_dim+1) # one extra feat: sequence length/ max_len
        x_bilstm = self.bilstm(seq_embedding)
        x_gearnet = self.gnn(graph_data, graph_data.residue_feature.float())
        print(x_bilstm.shape, x_gearnet.shape)

        x_integrated = torch.cat((x_bilstm, x_gearnet), dim=1).unsqueeze(1).transpose(1, 2) 

        scores = self.attention_fc(x_integrated) # attention across feature dimensions
        alpha = F.softmax(scores, dim=1)
        x_integrated = torch.sum(alpha * x_integrated, dim=1)

        x = F.relu(self.fc1(x_integrated)) # fully connected layer
        x = F.sigmoid(self.fc2(x))
        return x.squeeze(-1)
    
class GCN(nn.Module):
    def __init__(self, args):
        super(GCN, self).__init__()
        self.feat_in_dim = args['feat_in_dim']
        self.topo_in_dim = args['topo_in_dim']
        self.dims = [self.feat_in_dim+self.topo_in_dim] + list(args['hidden_dims'])
        self.module_out_dim = args['GNN_out_dim']
        self.batch_size = args['batch_size']

        self.layers = nn.ModuleList()
        for i in range(len(self.dims) - 1):
            self.layers.append(layers.GraphConv(self.dims[i], self.dims[i + 1], activation='relu'))

        self.batch_norms = nn.ModuleList()
        for i in range(len(self.dims) - 1):
            self.batch_norms.append(nn.BatchNorm1d(self.dims[i + 1]))

        self.mlp = layers.MLP(args['hidden_dims'][-1], 1,
            batch_norm=False, dropout=0, activation = 'relu')

        self.readout = layers.SumReadout()
        self.dropout = nn.Dropout(p=0.5)

    def forward(self, graph, input):
        hiddens = []
        layer_input = input

        for i in range(len(self.layers)):
            hidden = self.layers[i](graph, layer_input)
            hidden = self.batch_norms[i](hidden)
            hidden = self.dropout(hidden)
            hiddens.append(hidden)
            layer_input = hidden
        
        node_feature = hiddens[-1]
        graph_feature = self.readout(graph, node_feature)
        pred = self.mlp(graph_feature).squeeze(-1)

        return nn.functional.sigmoid(pred)