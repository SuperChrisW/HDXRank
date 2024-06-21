import os
import numpy as np
import pandas as pd
from collections import defaultdict
from pepGraph_dataset import pepGraph, pep

from torch_geometric.data import Dataset, Data
import torch
from tqdm import tqdm

cluster_index = 2

if __name__ == '__main__':
    root_dir = '/home/lwang/models/HDX_LSTM/data/Latest_test'
    save_dir = f'{root_dir}/graph_ensemble_GearNetEdge/cluster{cluster_index}_hop0'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    df = pd.read_excel(f'{root_dir}/merged_data.xlsx', sheet_name='Sheet2')
    df = df.dropna(subset=['structure_file'])
    apo_identifier = df['structure_file'].astype(str).unique()

    protein, state, chain_identifier = [], [], []
    correction = []
    match_uni = []
    database_id = []
    protein_chains=  []
    complex_state = []

    #process HDX data by apo_identifier (pdb structures)
    for i, temp_apo in enumerate(apo_identifier):
        temp_df = df[(df['structure_file'] == temp_apo)]
        temp_protein = temp_df['protein'].astype(str).to_list()
        temp_state = temp_df['state'].astype(str).to_list()
        temp_chain = temp_df['chain_identifier'].astype(str).to_list()
        temp_correction = temp_df['correction_value'].astype(int).to_list()
        temp_uni = temp_df['match_uni'].astype(str).to_list()
        temp_protein_chains= temp_df['protein_chain'].astype(str).to_list()
        temp_complex_state = temp_df['complex_state'].astype(str).to_list()

        protein.append(temp_protein)
        state.append(temp_state)
        chain_identifier.append(temp_chain)
        correction.append(temp_correction)
        match_uni.append(temp_uni)
        database_id.extend(temp_df['database_id'].astype(str).unique())
        protein_chains.append(temp_protein_chains[0])
        complex_state.append(temp_complex_state[0])


    keys = [database_id, protein, state, match_uni, apo_identifier, chain_identifier, correction, protein_chains, complex_state]
    print('total number of keys:', len(database_id))
    print(len(apo_identifier))

    graph_dataset = pepGraph(keys, root_dir, cluster_index, nfeature = 56, min_distance = 5.0, max_distance = 10.0,
                             truncation_window_size = None) # distance cutoff: distance between CA atoms, consider 10 or 15
    
    count = 0
    progress = tqdm(enumerate(graph_dataset), total=len(graph_dataset))
    for i, data in progress:
        if data is None:
            continue
        graph_ensemble, label = data
        path = f'{save_dir}/{label}.pt'
        torch.save(graph_ensemble, path)
        count += len(graph_ensemble)

    print(count)


