import argparse
import pickle
import os
import numpy as np
import pandas as pd
from collections import defaultdict
from pepGraph_dataset import pepGraph, pep

from torch_geometric.data import Dataset, Data
import torch
from tqdm import tqdm

if __name__ == '__main__':
    root_dir = '/home/lwang/models/HDX_LSTM/data/COVID_SPIKE'
    save_dir = f'{root_dir}/graph_ensemble/v2_ensemble'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    df = pd.read_excel(f'{root_dir}/COVID_record.xlsx', sheet_name='Sheet1')
    df = df.dropna(subset=['chain_identifier'])
    apo_identifier = df['apo_identifier'].astype(str).unique()

    protein, state, chain_identifier = [], [], []
    correction = []
    match_uni = []
    database_id = []
    protein_chains=  []

    #process HDX data by apo_identifier (pdb structures)
    for temp_apo in apo_identifier:
        temp_df = df[(df['apo_identifier'] == temp_apo)]
        temp_protein = temp_df['protein'].astype(str).to_list()
        temp_state = temp_df['state'].astype(str).to_list()
        temp_chain = temp_df['chain_identifier'].astype(str).to_list()
        temp_correction = temp_df['correction_value'].astype(int).to_list()
        temp_uni = temp_df['match_uni'].astype(str).to_list()
        temp_protein_chains= temp_df['protein_chain'].astype(str).to_list()
        print(temp_apo, temp_protein_chains)

        protein.append(temp_protein)
        state.append(temp_state)
        chain_identifier.append(temp_chain)
        correction.append(temp_correction)
        match_uni.append(temp_uni)
        database_id.extend(temp_df['database_id'].astype(str).unique())
        protein_chains.append(temp_protein_chains[0])
        
    keys = [database_id, protein, state, match_uni, apo_identifier, chain_identifier, correction, protein_chains]
    print('total number of keys:', len(database_id))
    print(len(apo_identifier))

    graph_dataset = pepGraph(keys, root_dir)
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


