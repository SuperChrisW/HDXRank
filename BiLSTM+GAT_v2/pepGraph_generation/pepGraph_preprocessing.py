import os
import numpy as np
import pandas as pd
from collections import defaultdict
from pepGraph_dataset import pepGraph, pep

from torch_geometric.data import Dataset, Data
import torch
from tqdm import tqdm

cluster_index = 2
protein_model = '8F7A'

if __name__ == '__main__':
    root_dir = '/home/lwang/models/HDX_LSTM/data/Latest_set'
    #save_dir = f'{root_dir}/graph_ensemble_GVP/{protein_model}/cluster{cluster_index}'
    save_dir = f'{root_dir}/graph_ensemble_modifiedGVP/cluster{cluster_index}'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    df = pd.read_excel(f'{root_dir}/merged_data.xlsx', sheet_name='Sheet1')
    df = df.dropna(subset=['structure_file'])
    #df = df[df['note'] == protein_model]

    apo_identifier = list(df['structure_file'].astype(str).unique())
    protein, state, chain_identifier = [], [], []
    correction = []
    match_uni = []
    database_id = []
    protein_chains=  []
    complex_state = []
    structure_list = []

    #process HDX data by apo_identifier (pdb structures)
    for i, temp_apo in enumerate(apo_identifier):
        if temp_apo.split(":")[0] != 'MODEL':
            structure_list.append(temp_apo)
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
        else:
            model_list = [f'MODEL_{i}_REVISED' for i in range(1, 501)]
            model = temp_apo.split(":")[1]
            temp_df = df[df['structure_file'] == model]

            temp_model_collection = pd.DataFrame()
            for temp_model in model_list:
                temp_df.loc[:, 'structure_file'] = temp_model

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

            structure_list = structure_list + model_list

    keys = [database_id, protein, state, match_uni, structure_list, chain_identifier, correction, protein_chains, complex_state]
    print('total number of keys:', len(database_id))

    graph_dataset = pepGraph(keys, root_dir, cluster_index, nfeature = 56, min_distance = 5.0, max_distance = 10.0, protein_name = protein_model,
                             truncation_window_size = None, graph_type='GVP') # graph_type: switch between 'GearNet' and 'GVP'
    
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


