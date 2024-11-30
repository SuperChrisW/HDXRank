import os
import numpy as np
import pandas as pd
from collections import defaultdict
from pepGraph_dataset import pepGraph, pep

from torch_geometric.data import Dataset, Data
import torch
from tqdm import tqdm

def graph_tasks(df, N_model=1):
    apo_identifier = list(df['structure_file'].astype(str).unique())
    protein, state, chain_identifier = [], [], []
    correction = []
    database_id = []
    protein_chains=  []
    complex_state = []
    structure_list = []
    feat_identifier = []

    #process HDX data by apo_identifier (pdb structures)
    for i, temp_apo in enumerate(apo_identifier):
        if temp_apo.split(":")[0] != 'MODEL':
            temp_df = df[(df['structure_file'] == temp_apo)]
            temp_protein = temp_df['protein'].astype(str).to_list()
            temp_state = temp_df['state'].astype(str).to_list()
            temp_chain = temp_df['chain_identifier'].astype(str).to_list()
            temp_correction = temp_df['correction_value'].astype(int).to_list()
            temp_protein_chains= temp_df['protein_chain'].astype(str).to_list()
            temp_complex_state = temp_df['complex_state'].astype(str).to_list()

            structure_list.append(temp_apo)
            feat_identifier.append([temp_apo.upper()]*len(temp_protein_chains))
            protein.append(temp_protein)
            state.append(temp_state)
            chain_identifier.append(temp_chain)
            correction.append(temp_correction)
            database_id.extend(temp_df['database_id'].astype(str).unique())
            protein_chains.append(temp_protein_chains)
            complex_state.append(temp_complex_state[0])
        else:
            model_list = [f'MODEL_{i}_REVISED' for i in range(1, N_model+1)]
            apo_models = temp_apo.split(":")[1:] # suppose the format is MODEL:apo1:apo2: ...
            temp_df = df[df['structure_file'].isin(apo_models)]

            temp_protein = [temp_df['protein'].astype(str).to_list()] * N_model
            temp_state = [temp_df['state'].astype(str).to_list()] * N_model
            temp_chain = [temp_df['chain_identifier'].astype(str).to_list()] * N_model
            temp_correction = [temp_df['correction_value'].astype(int).to_list()] * N_model
            temp_complex_state = ['protein complex'] * N_model
            temp_database_id = [temp_df['database_id'].astype(str).to_list()[0]] * N_model

            temp_protein_chains= [temp_df['protein_chain'].astype(str).to_list()] * N_model

            protein.extend(temp_protein)
            state.extend(temp_state)
            chain_identifier.extend(temp_chain)
            correction.extend(temp_correction)
            database_id.extend(temp_database_id)
            protein_chains.extend(temp_protein_chains)
            complex_state.extend(temp_complex_state)
            feat_identifier.extend([apo_models] * N_model)

            structure_list = structure_list + model_list
    print('total number of keys:', len(database_id))
    keys = [database_id, protein, state, structure_list, chain_identifier, correction, protein_chains, complex_state, feat_identifier]
    return keys

def main(task_table, usage, graph_type, embedding_type, radius_max, radius_min, cluster_index, protein_model, N_model):
    if usage == 'train':
        root_dir = '/home/lwang/models/HDX_LSTM/data/Latest_test'
        save_dir = f'{root_dir}/graph_ensemble_simple{graph_type}/cluster{cluster_index}_{radius_max}A_{embedding_type}_rescale/'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        #df = pd.read_excel(f'{root_dir}/merged_data-NAcomplex.xlsx', sheet_name='Sheet1')
        df = pd.read_excel(f'{root_dir}/test_data.xlsx', sheet_name='Sheet2')
        df = df.dropna(subset=['structure_file'])
        df = df[df['complex_state'] != 'ligand complex']
    elif usage == 'hdock':
        root_dir = '/home/lwang/models/HDX_LSTM/data/hdock'
        save_dir = f'{root_dir}/graph_ensemble_{graph_type}/{protein_model}/cluster{cluster_index}_{radius_max}A_manual_rescale'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        df = pd.read_excel(f'{root_dir}/test_data_AF.xlsx', sheet_name='hdock') ## 'hdock_ori' or 'hdock_modeling'
        df = df.dropna(subset=['structure_file'])
        df = df[df['note'] == protein_model[:4]]

    elif usage == 'covid':
        root_dir = '/home/lwang/models/HDX_LSTM/data/COVID_SPIKE'
        save_dir = f'{root_dir}/graph_ensemble_{graph_type}/{protein_model}/cluster{cluster_index}_{radius_max}A_manual_rescale'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        df = pd.read_excel(f'{root_dir}/{task_table}.xlsx', sheet_name='Sheet1')
        df = df.dropna(subset=['structure_file'])
        df = df[df['note'] == protein_model]

    tasks = graph_tasks(df, N_model=N_model)
    graph_dataset = pepGraph(tasks, root_dir, save_dir, cluster_index, min_distance = radius_min, max_distance = radius_max, protein_name = protein_model,
                            graph_type=graph_type, embedding_type=embedding_type, usage=usage) 

    count = 0
    progress = tqdm(enumerate(graph_dataset), total=len(graph_dataset))
    for i, data in progress:
        if data is None:
            continue
        graph_ensemble, label = data
        path = f'{save_dir}/{label}.pt'
        if len(graph_ensemble) == 0:
            continue
        torch.save(graph_ensemble, path)
        count += len(graph_ensemble)
    print(count)

if __name__ == '__main__':
    cluster_index = 1
    usage = 'covid' # 'train' or 'hdock'
    graph_type = 'GearNet' # 'GearNet' or 'GVP'
    embedding_type = 'manual' # 'esm2' or 'manual'
    radius_max = 8
    radius_min = 3

    N_model = 1000
    spike_protein = 'Delta'
    proteins = [f'fold_vh16_vl106_seed2_model_0',
                f'fold_vh16_vl106_seed4_model_0',
                f'fold_vh16_vl106_seed4_model_1',
                f'fold_vh16_vl106_seed4_model_2',
                f'fold_vh16_vl106_seed5_model_2']
    proteins = [f'{spike_protein}_{protein}' for protein in proteins]
    task_table = f'COVID_{spike_protein}'

    for protein in proteins[4:5]:
        main(task_table, usage, graph_type, embedding_type, radius_max, radius_min, cluster_index, protein, N_model)

