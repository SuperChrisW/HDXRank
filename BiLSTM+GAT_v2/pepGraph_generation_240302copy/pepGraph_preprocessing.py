import argparse
import pickle
import os
import numpy as np
import pandas as pd
from collections import defaultdict
from pepGraph_dataset import pepGraph, pep

from torch_geometric.data import Dataset, Data
import torch

nfeatures = 48
def data_process(pep_registry, matrix):
    edge_index = []
    edge_attr = []
    n = matrix.shape[0]
    embedding = torch.zeros((n, nfeatures), dtype=torch.float32) # may include multiple chains, n = total(peptides)
    range_mtx = torch.zeros((n, 2), dtype=torch.int)
    hdx_value = torch.zeros((n, 1), dtype=torch.float32)
    index_dict = {}
    pdist = torch.nn.PairwiseDistance(p=2)

    for i in range(n):
        edges = [[i, j] for j, val in enumerate(matrix[i]) if val == 1]
        xyz1 = torch.tensor([pep_registry[i].position] * len(edges))
        xyz2 = torch.tensor([pep_registry[j].position for j, val in enumerate(matrix[i]) if val == 1])
        d2 = pdist(xyz1, xyz2).reshape(-1)

        edge_index.extend(edges)
        #edge_attr.extend([1 / val if val > 1 else 1 for val in d2])
        edge_attr.extend(d2)
    edge_attr = torch.tensor(edge_attr, dtype=torch.float32)

    for i, pep in enumerate(pep_registry):
        embedding[i] = pep.embedding
        range_mtx[i] = torch.tensor([pep.start, pep.end])
        index_dict[f'{pep.chain}_{pep.start}_{pep.end}'] = pep.i
        hdx_value = pep.hdx_value

    return embedding, torch.tensor(edge_index, dtype=torch.int), range_mtx, index_dict, hdx_value, edge_attr

root_dir = '/home/lwang/AI-HDX-main/HDX_MS_dataset/complexpair_dataset'
df = pd.read_excel(f'{root_dir}/merged_complex_pair.xlsx', sheet_name='Sheet1')
df = df.dropna(subset=['chain_identifier'])
apo_identifier = df['apo_identifier'].astype(str).unique()

protein, state, chain_identifier = [], [], []
correction = []
match_uni = []
database_id = []
for temp_apo in apo_identifier:
    temp_df = df[(df['apo_identifier'] == temp_apo)]
    temp_protein = temp_df['protein'].astype(str).to_list()
    temp_state = temp_df['state'].astype(str).to_list()
    temp_chain = temp_df['chain_identifier'].astype(str).to_list()
    temp_correction = temp_df['correction_value'].astype(int).to_list()
    temp_uni = temp_df['match_uni'].astype(str).to_list()
    protein.append(temp_protein)
    state.append(temp_state)
    chain_identifier.append(temp_chain)
    correction.append(temp_correction)
    match_uni.append(temp_uni)
    database_id.extend(temp_df['database_id'].astype(str).unique())


#database_id = df['database_id'].astype(str).to_list()
#apo_identifier = df['apo_identifier'].astype(str).to_list()
keys = [database_id, protein, state, match_uni, apo_identifier, chain_identifier, correction]
print('total number of keys:', len(database_id))
print(len(apo_identifier))

train_dataset = pepGraph(keys, root_dir)

count = 0

for i, data in enumerate(train_dataset):
    pep_graph = None
    if data is None:
        continue
    label, pep_registry, pep_adj = data

    x, edge_index, range_mtx, index_dict, hdx_value, edge_attr = data_process(pep_registry, pep_adj)
    path = f'{root_dir}/pepGraph_files/{label}.pt'

    #if os.path.exists(path):
    #    continue
    pep_graph = Data(x=x, y=hdx_value, edge_index = edge_index.t().contiguous())
    pep_graph['range'] = range_mtx
    pep_graph['pep_index'] = index_dict
    pep_graph['embedding'] = pep.protein_embedding # dict[chain] = embedding_matrix (N(res), 46)
    pep_graph.edge_attr = edge_attr

    pep_graph.validate(raise_on_error=True)
    torch.save(pep_graph, path)
    count += len(pep_registry)

print(count)


