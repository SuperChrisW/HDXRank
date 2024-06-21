from torchdrug import data, utils, tasks, core, layers
from torchdrug.layers import geometry
from GearNet import GearNet
from itertools import islice
import torch

pdb_fpath = '/home/lwang/models/HDX_LSTM/data/Fullset/structure/AF-O60664-F1-model_v4.pdb'
protein = data.Protein.from_pdb(pdb_fpath)
protein.view = 'residue'
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

graph_construction_model = layers.GraphConstruction(node_layers=[geometry.AlphaCarbonNode()], 
                                                    edge_layers=[geometry.SpatialEdge(radius=10.0, min_distance=5),
                                                                 geometry.KNNEdge(k=10, min_distance=5),
                                                                 geometry.SequentialEdge(max_distance=2)])

_protein = data.Protein.pack([protein]*4)
protein_ = graph_construction_model(_protein)
print("Graph before: ", _protein)
print("Graph after: ", protein_)

model = GearNet(input_dim = 67, hidden_dims = [512],
                    num_relation=7, batch_norm=True, concat_hidden=True, readout='sum', short_cut=True).to(device)


train_loader = data.DataLoader(protein_, batch_size=4)
for i, batch in enumerate(islice(train_loader, len(train_loader))):
    batch = batch.to(device)
    print(batch)
    output = model(batch, batch.node_feature.float())
    print(output.shape)
