"""
2024/3/28
Author: WANG Liyao
Note: 
Take the HHBlits output hhm format, PDB structure as input to encode peptide sequence and structure information into embedding.
The Heteroatoms (NA and SM) encoding are also supported and merged with protein sequence,
but it is not used in the current HDXRank model.
"""
import os
import torch
import torch.nn.functional as F
import numpy as np
import pandas as pd
from Bio import PDB
from Bio.PDB import NeighborSearch, Selection
import warnings
from tool.BioWrappers import get_bio_model
from pepGraph_utlis import load_protein, load_nucleic_acid, load_sm, RawInputData, ChemData

def find_contact_res(HETATM_input, res_list, cutoff):
    '''
    find the contact residues of the heteroatoms(NA or SM)
    return contact matrix: [#res, #entity of NA/SM] where elements are type encoding
    '''
    entity_index_map = {tuple(entity['coord']): entity['token_type'] for entity in HETATM_input.seq_data.values()}
    res_index_map = {res: index for index, res in enumerate(res_list)}
    contact_mtx = np.full((len(res_list), len(entity_index_map.keys())), 20) # 20 for 'UNK' unknown type

    atom_list = Selection.unfold_entities(res_list, 'A')
    ns = NeighborSearch(atom_list)
    for entity_idx, entity_coord in enumerate(entity_index_map.keys()):
        for nearby_res in ns.search(entity_coord, cutoff, level='R'): # level='R' means search by residue
            res_idx = res_index_map[nearby_res]
            contact_mtx[res_idx, entity_idx] = entity_index_map[entity_coord] # fill entity type into contact matrix

    return contact_mtx

def merge_inputs(inputs_list):
    """
    Merge the input data from different chains
    """
    if len(inputs_list) == 0:
        return RawInputData()
    elif len(inputs_list) == 1:
        return inputs_list[0]
    else:
        running_input = inputs_list[0] # RawInputData object
        for i in range(1, len(inputs_list)):
            running_input = running_input.merge(inputs_list[i])
        return running_input

def embed_protein(structure_dir, hhm_dir, save_dir, pdb_fname, protein_chain_hhms, NA_chain = [], SM_chain = []):
    print('processing:', pdb_fname)
    protein_chain = list(protein_chain_hhms.keys())
    if len(protein_chain) == 0:
        raise ValueError('protein_chain is empty')

    pdb_file = f'{structure_dir}/{pdb_fname}.pdb'
    if not os.path.isfile(pdb_file):
        print('file not exist:', pdb_file)
        return None
    structure = get_bio_model(pdb_file)
    chains = list(structure.get_chains())
    print('protein chains:', protein_chain)
    print('NA chains:', NA_chain)
    print('SM chains:', SM_chain)

    ### load input data from protein and heteroatoms ###
    protein_inputs = []
    NA_inputs = []
    SM_inputs = []

    for chain in chains:
        chain_id = chain.get_id()
        hhm_file = os.path.join(hhm_dir, f'{protein_chain_hhms[chain_id]}_{chain_id}.hhm')
        if chain_id in protein_chain:
            if os.path.isfile(hhm_file) and os.path.isfile(pdb_file):
                protein_inputs.append(load_protein(hhm_file, pdb_file, chain_id))
            else:
                raise ValueError('file not exist', hhm_file, pdb_file, chain_id)
        elif chain_id in NA_chain:
            NA_inputs.append(load_nucleic_acid(pdb_file, chain_id))
        elif chain_id in SM_chain:
            SM_inputs.append(load_sm(pdb_file, chain_id))

    ### merging protein, NA, and SM ###
    embedding = [protein.construct_embedding() for protein in protein_inputs]
    embed_mtx = torch.cat(embedding, dim=0)

    ### find contact residues of NA and SM ###
    chain_list = [chain for chain in structure.get_chains() if chain.id in protein_chain]
    res_list = Selection.unfold_entities(chain_list, 'R')

    ### block the following code as the Heteroatom encoding is not used in the current HDXRank model ###
    '''res_idx_list = []
    for res in res_list:
        name = res.get_resname()
        res_idx_list.append(chemdata.aa2num[name] if name in chemdata.aa2num else 20)
    res_idx_list = np.array(res_idx_list).reshape(-1, 1)

    contact_ensemble = []
    if len(NA_inputs) != 0 or len(SM_inputs) != 0:
        for inputs in [NA_inputs, SM_inputs]:
            merged_HETinput = merge_inputs(inputs)
            if len(merged_HETinput.seq_data.keys()) == 0:
                continue
            contact_mtx = find_contact_res(merged_HETinput, res_list, cutoff = 5.0)  # [#res, #entity of NA/SM] where elements are type encoding
            contact_ensemble.append(contact_mtx)

    contact_ensemble.insert(0, res_idx_list)
    contact_mtx = np.concatenate(contact_ensemble, axis=1) if len(contact_ensemble) > 1 else contact_ensemble[0]

    contact_tensor = torch.tensor(contact_mtx, dtype=torch.long).flatten()
    encoded_tensor = F.one_hot(contact_tensor, num_classes=len(chemdata.num2aa))

    encoded_tensor = encoded_tensor.view(contact_mtx.shape[0], -1, encoded_tensor.shape[1])
    encoded_tensor = torch.sum(encoded_tensor, dim=1)
    encoded_tensor = torch.log(encoded_tensor + 1) # apply log(e) to elements

    protein_embedding = torch.cat((embed_mtx, encoded_tensor), dim=1)
    print('protein_embedding:', protein_embedding.shape)'''

    #save to file
    res_idx = [res.id[1] for res in res_list]
    res_name = [res.get_resname() for res in res_list]
    chain_label = [res.get_parent().id for res in res_list]
    data_to_save = {
        'res_idx': res_idx,
        'res_name': res_name,
        'chain_label': chain_label,
        'embedding': embed_mtx
    }
    save_file = os.path.join(save_dir, f'{pdb_fname}.pt')
    torch.save(data_to_save, save_file)
    
def generate_embedding(root_dir, task_table, protein_name=None, N_model=None):
    warnings.filterwarnings("ignore")
    hhm_dir = os.path.join(root_dir, 'hhm_files', protein_name)
    save_dir = os.path.join(root_dir, 'embedding_files', protein_name)
    structure_dir = os.path.join(root_dir, 'structure', protein_name)
    df = pd.read_excel(f'{root_dir}/{task_table}.xlsx', sheet_name='Sheet1')
    df = df[df['note'] == protein_name] if protein_name is not None else df # filter, embedding only for the specified protein
    df = df.dropna(subset=['structure_file']).drop_duplicates(subset=['structure_file'])

    for id, row in df.iterrows():
        pdb_fname = str(row['structure_file']).upper().split('.')[0]

        protein_chain = row['protein_chain'].split(',')
        NA_chain = row['NA_chain'].split(',') if not pd.isnull(row['NA_chain']) else []
        SM_chain = row['SM_chain'].split(',') if not pd.isnull(row['SM_chain']) else []

        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        if os.path.isfile(f'{save_dir}/{pdb_fname}.pt'):
            print('file exist')
            continue

        protein_chain_hhms = {}
        if pdb_fname.split(":")[0] != 'MODEL':
            if row['complex_state'] == 'protein complex':
                continue
            for i, chain in enumerate(protein_chain):
                protein_chain_hhms[chain] = pdb_fname
            embed_protein(structure_dir, hhm_dir, save_dir, pdb_fname,
                                protein_chain_hhms = protein_chain_hhms, NA_chain = NA_chain, SM_chain = SM_chain)
        else:
            #hhm_fname should be a list of hhm file names, corresponding to protein_chain
            model_list = [f'MODEL_{i}_REVISED' for i in range(1, N_model+1)]
            hhm_fname = row['structure_file'].split(':')[1:]
            if len(hhm_fname) != len(protein_chain):
                raise ValueError('hhm_fname and protein_chain do not match')
            for i, chain in enumerate(protein_chain):
                protein_chain_hhms[chain] = hhm_fname[i].upper()

            for model in model_list:
                embed_protein(structure_dir, hhm_dir, save_dir, model,
                                    protein_chain_hhms = protein_chain_hhms, NA_chain = NA_chain, SM_chain = SM_chain)

if __name__ == "__main__":

    ## covid spike protein
    '''root_dir = '/home/lwang/models/HDX_LSTM/data/COVID_SPIKE/'
    spike_protein = 'Delta'
    proteins = [f'fold_vh16_vl106_seed2_model_0',
                f'fold_vh16_vl106_seed4_model_0',
                f'fold_vh16_vl106_seed4_model_1',
                f'fold_vh16_vl106_seed4_model_2',
                f'fold_vh16_vl106_seed5_model_2']
    proteins = [f'{spike_protein}_{protein}' for protein in proteins]'''

    # training data
    '''root_dir = '/home/lwang/models/HDX_LSTM/data/Latest_set'
    hhm_dir = os.path.join(root_dir, 'hhm_files')
    save_dir = os.path.join(root_dir, 'embedding_files')
    structure_dir = os.path.join(root_dir, 'structure')
    df = pd.read_excel(f'{root_dir}/merged_data-NAcomplex.xlsx', sheet_name='Sheet1')
    df = df.dropna(subset=['structure_file']).drop_duplicates(subset=['structure_file'])'''