import torch
from torchdrug import data

from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error
from math import sqrt
import os
import sys

sys.path.append("/home/lwang/AI-HDX-main/ProteinComplex_HDX_prediction/BiLSTM+GAT_v2/pepGraph_generation")
from GearNet import GearNet

def test_model(model, test_loader, device, rm_feat):
    y_pred = []
    y_true = []
    model.eval()
    with torch.no_grad():
        for graph_batch in test_loader: 
            graph_batch = graph_batch.to(device)
            targets = graph_batch.y
            if rm_feat == 'msa':
                node_feat = graph_batch.residue_feature[:,30:].float()
            elif rm_feat == 'sequence':
                node_feat = torch.cat([graph_batch.residue_feature[:,:30].float(), graph_batch.residue_feature[:,40:].float()], dim=1) # remove seq feat
            elif rm_feat == 'physical':
                node_feat = torch.cat([graph_batch.residue_feature[:,:40].float(), graph_batch.residue_feature[:,44:].float()], dim=1) # remove physical feat
            elif rm_feat == 'geometric':
                node_feat = torch.cat([graph_batch.residue_feature[:,:44].float(), graph_batch.residue_feature[:,56:].float()], dim=1) # remove geometric feat
            elif rm_feat == 'heteroatom':
                node_feat = graph_batch.residue_feature[:,:56].float() # remove heteroatom feat
            elif rm_feat == 'none':
                node_feat = graph_batch.residue_feature.float()
            elif rm_feat == 'esm2':
                node_feat = graph_batch.residue_feature.float()
            #node_feat = torch.cat([graph_batch.residue_feature[:,:35].float(), graph_batch.residue_feature[:,40:41].float()], dim=1) # same as AI-HDX
            outputs = model(graph_batch, node_feat)

            y_pred.extend(outputs.cpu().detach().numpy())
            y_true.extend(targets.cpu().detach().numpy())

    return y_true, y_pred

##################################### initial setting ##################################### 
root_dir = '/home/lwang/models/HDX_LSTM/data/Latest_test'
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
cluster = 'cluster1_8A_esm2_rescale'
rm_feat = 'esm2'
model_name = 'GearNet'
model_fpath = f'/home/lwang/models/HDX_LSTM/results/241110_GearNet/model_GN56_{cluster}_v0_{rm_feat}_epoch99.pth'
pepGraph_dir = os.path.join(root_dir, 'graph_ensemble_simpleGearNet', cluster)

##################################### data loading ##################################### 
input_graph = []
for file in os.listdir(pepGraph_dir):
    pepGraph_file = f'{pepGraph_dir}/{file}'
    pepGraph_ensemble = torch.load(pepGraph_file)
    input_graph.extend(pepGraph_ensemble)
print(f"Total number of graphs: {len(input_graph)}")

#GearNet
feat_num = {"sequence": 10, "msa": 30, "physical": 4, "geometric": 12, "heteroatom": 42, 'none': 0, '36feat': 62, 'esm2': -264}
model = GearNet(input_dim = 56-feat_num[rm_feat], hidden_dims = [512,512,512],
                num_relation=7, batch_norm=True, concat_hidden=True, readout='sum', activation = 'relu', short_cut=True)

model_state_dict = torch.load(model_fpath, map_location=device)
model_state_dict = model_state_dict['model_state_dict'] # add if saved as checkpoint
model.load_state_dict(model_state_dict)
model = model.to(device)

total_set = data.Protein.pack(input_graph)
total_dataloader = data.DataLoader(total_set, batch_size=16, shuffle=False)

##################################### evaluation ##################################### 
y_true, y_pred = test_model(model, total_dataloader, device, rm_feat)
rmse = sqrt(mean_squared_error(y_true, y_pred))
pcc = pearsonr(y_true, y_pred)[0]
spR = spearmanr(y_true, y_pred)[0]
print(f'PCC: {pcc}, SPR: {spR}, RMSE: {rmse}')
