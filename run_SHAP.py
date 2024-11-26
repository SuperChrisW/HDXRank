import numpy as np
from itertools import chain, combinations
from math import factorial
import random

import torch
import torchdrug as td
import os
import sys
sys.path.append('/home/lwang/AI-HDX-main/ProteinComplex_HDX_prediction/HDXRank')
from tqdm import tqdm
from GearNet import GearNet

def mask_features_with_mean(X, fixed_features, baseline = 'graph', GLOBAL_FEATURE_MEANS = None):
    """
    Replace unfixed features in the input feature matrix with their mean values across the dataset.

    Parameters:
    - X: np.ndarray, shape (N_nodes, N_dims)
        The input feature matrix.
    - fixed_features: list of int
        List of feature indices (in the last dimension) to be kept as fixed.
    
    Returns:
    - masked_X: np.ndarray, shape (N_nodes, N_dims)
        The feature matrix with unfixed features replaced by their respective mean values.
    """
    # Validate input
    if not isinstance(X, (np.ndarray, torch.Tensor)):
        raise ValueError("Input X must be a numpy array or torch tensor.")
    if any(idx >= X.shape[-1] or idx < 0 for idx in fixed_features):
        raise ValueError("Invalid feature indices in fixed_features.")

    # Compute the mean values of each feature across the dataset
    if baseline == 'graph':
        feature_means = np.mean(X, axis=0)  # Shape: (, N_dims)
    elif baseline == 'dataset':
        feature_means = GLOBAL_FEATURE_MEANS

    # Replace unfixed features with their mean values
    masked_X = np.copy(X)
    for feature_idx in range(X.shape[1]):
        if feature_idx not in fixed_features:
            masked_X[:, feature_idx] = feature_means[feature_idx]

    return masked_X

def generate_coalitions(feat_set):
    """
    Generate all possible coalitions (subsets) of fixed features for a given number of features.

    Parameters:
    - N_wo_j: int
        total feature index set (N) except for the feature index j.

    Returns:
    - coalitions: list of lists
        A nested list containing all subsets of feature indices, including the empty set and the full set.
    """
    # Generate all possible subsets (power set)
    N_len = len(feat_set)
    coalitions = list(chain.from_iterable(combinations(feat_set, r) for r in range(N_len + 1)))
    
    # Convert each subset (tuple) to a list for compatibility with mask_features_with_mean
    coalitions = [list(coalition) for coalition in coalitions]
    
    return coalitions

def compute_coalition_prob(S_len, N):
    """
    Compute the probability of each coalition size for a given number of features.

    Parameters:
    - S_len: int
        number of fixed features

    Returns:
    """
    prob = (factorial(S_len) * factorial(N - S_len - 1)) / factorial(N)
    return prob

def shapley_value_GNN(graph, node_embedding, target_feat, model, GNN_feats, GLOBAL_FEATURE_MEANS):
    graph = graph.to(device)
    model = model.to(device)

    N_features = node_embedding[0].shape[-1]
    N_wo_j = list(GNN_feats.keys())
    N_wo_j.remove(target_feat)

    S_coalitions = generate_coalitions(N_wo_j)
    S_probs = {S_len: compute_coalition_prob(S_len, N_features) for S_len in range(N_features)}

    shapley = 0
    model.eval()
    for S in S_coalitions:
        # convert S to real feature index
        real_S = []
        for feat in S:
            real_S.extend(GNN_feats[feat])
        real_S_plusJ = real_S + GNN_feats[target_feat]

        masked_X1 = mask_features_with_mean(node_embedding, real_S, baseline='dataset', GLOBAL_FEATURE_MEANS=GLOBAL_FEATURE_MEANS)
        masked_X2 = mask_features_with_mean(node_embedding, real_S_plusJ, baseline='dataset', GLOBAL_FEATURE_MEANS=GLOBAL_FEATURE_MEANS)
        masked_X1_tensor = torch.tensor(masked_X1).to(device)
        masked_X2_tensor = torch.tensor(masked_X2).to(device)

        marginal_contribution = model(graph, masked_X2_tensor).cpu().detach().numpy() - model(graph, masked_X1_tensor).cpu().detach().numpy()
        shapley += marginal_contribution * S_probs[len(S)]
    return shapley.item()


## GNN graphs loading
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

data_dir = '/home/lwang/models/HDX_LSTM/data/Latest_set/graph_ensemble_simpleGearNet/cluster1_8A_manual_rescale'
graphs = []
for i, file in tqdm(enumerate(os.listdir(data_dir))):
    if file.endswith('.pt'):
        graphs.extend(torch.load(os.path.join(data_dir, file)))
print(len(graphs))

## GNN model loading
model_fpath = f'/home/lwang/models/HDX_LSTM/results/240918_GVP/model_GN56_cluster1_8A_manual_rescale_v0_epoch99.pth'
model = GearNet(56, hidden_dims = [512,512,512], num_relation=7, batch_norm=True, concat_hidden=True, readout='sum', activation = 'relu', short_cut=True)
model_state_dict = torch.load(model_fpath, map_location=device)
model_state_dict = model_state_dict['model_state_dict'] # add if saved as checkpoint
model.load_state_dict(model_state_dict)

## Compute SHAP values
random.seed(42)
sample_graphs = random.sample(graphs, 2000)
sample_node_embedding = [graph.residue_feature.float() for graph in sample_graphs]
GLOBAL_FEATURE_MEANS = torch.mean(torch.cat(sample_node_embedding), dim=0)

GNN_feats = {
            "MSA": list(range(0, 30)),
            "HDMD": list(range(30, 35)),
            "Residue polarity": list(range(35, 39)),
            "Residue charge": list(range(39, 40)),
            "SASA": list(range(40, 41)),
            "HSE": list(range(41, 44)),
            "Dihedrals": list(range(44, 50)),
            "Orientations": list(range(50, 56))
            }


explanations = np.zeros((len(sample_graphs), len(GNN_feats.keys())))
for i in range(len(sample_graphs)):
    for j, feat_exp in tqdm(enumerate(GNN_feats.keys())):
        explanations[i, j] = shapley_value_GNN(sample_graphs[i], sample_node_embedding[i], feat_exp, model, GNN_feats, GLOBAL_FEATURE_MEANS)

np.save('/home/lwang/models/HDX_LSTM/results/241110_GearNet/SHAP_values_GearNet_cluster1_8A_manual_rescale_v0_epoch99.npy', explanations)