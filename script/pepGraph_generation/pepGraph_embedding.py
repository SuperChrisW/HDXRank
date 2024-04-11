"""
2024/3/28
Author: Liyao
Usage: Take the HHBlits output hhm format, and encode the protein sequence
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
    
if __name__ == "__main__":
    ## generate embedding ##
    warnings.filterwarnings("ignore")
    root_dir = '/home/lwang/models/HDX_LSTM/data/Fullset'

    dssp_dir = os.path.join(root_dir, 'dssp_files')
    hhm_dir = os.path.join(root_dir, 'hhm_files')
    rigidity_dir = os.path.join(root_dir, 'rigidity_files')
    save_dir = os.path.join(root_dir, 'embedding_files')
    structure_dir = os.path.join(root_dir, 'structrue')
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)

    df = pd.read_excel(f'{root_dir}/merged_data_oldVer.xlsx', sheet_name='Sheet1')
    df = df.dropna(subset=['chain_identifier'])
    df = df.drop_duplicates(subset=['apo_identifier'])
    fail_list = []
    chemdata = ChemData()

    for id, row in df.iterrows():
        apo_pdb = row['apo_identifier'].split('.')[0]
        print('processing:', apo_pdb)
        if os.path.isfile(f'{save_dir}/{apo_pdb}.pt'):
            print('file exist')
            continue
        if not pd.isnull(row['protein_chain']):
            protein_chain = row['protein_chain'].split(',')
        else:
            raise ValueError('protein_chain is empty')
        NA_chain = row['NA_chain'].split(',') if not pd.isnull(row['NA_chain']) else []
        SM_chain = row['SM_chain'].split(',') if not pd.isnull(row['SM_chain']) else []

        parser = PDB.PDBParser()
        pdb_file = f'{root_dir}/structure/{apo_pdb}.pdb'
        if not os.path.isfile(pdb_file):
            print('skip')
            continue
        structure = get_bio_model(pdb_file)
        chains = list(structure.get_chains())

        ### load input data from protein and heteroatoms ###
        protein_inputs = []
        NA_inputs = []
        SM_inputs = []

        try:
            for chain in chains:
                chain_id = chain.get_id()
                hhm_file = os.path.join(hhm_dir, f'{apo_pdb}_{chain_id}.hhm')
                rigidity_file = os.path.join(rigidity_dir, f'rigidity_{apo_pdb}.txt')

                if chain_id in protein_chain:
                    if os.path.isfile(hhm_file) and os.path.isfile(pdb_file):
                        protein_inputs.append(load_protein(hhm_file, rigidity_file, pdb_file, chain_id))
                    else:
                        fail_list.append((apo_pdb, chain_id))
                        raise ValueError('file not exist', apo_pdb, chain_id)
                elif chain_id in NA_chain:
                    NA_inputs.append(load_nucleic_acid(pdb_file, chain_id))
                elif chain_id in SM_chain:
                    SM_inputs.append(load_sm(pdb_file, chain_id))
        
            ### merging protein, NA, and SM ###
            embedding = [protein.construct_embedding() for protein in protein_inputs] # RawInputData object method
            embed_mtx = torch.cat(embedding, dim=0) #FIXME: not include index and residue type
        except Exception as e:
            print(e)
            continue

        ### find contact residues of NA and SM ###
        chain_list = [chain for chain in structure.get_chains() if chain.id in protein_chain]
        res_list = Selection.unfold_entities(chain_list, 'R') #contains all chains from protein

        res_idx_list = []
        for res in res_list:
            name = res.get_resname()
            res_idx_list.append(chemdata.aa2num[name] if name in chemdata.aa2num else 20)
        res_idx_list = np.array(res_idx_list).reshape(-1, 1)

        contact_ensemble = []        
        if len(NA_inputs) != 0 or len(SM_inputs) != 0:
            for inputs in [NA_inputs, SM_inputs]:
                merged_HETinput = merge_inputs(inputs)
                print('merged_HETinput:', merged_HETinput, type(merged_HETinput))
                if len(merged_HETinput.seq_data.keys()) == 0:
                    continue
                contact_mtx = find_contact_res(merged_HETinput, res_list, cutoff = 5.0)  # [#res, #entity of NA/SM] where elements are type encoding
                contact_ensemble.append(contact_mtx)

        contact_ensemble.insert(0, res_idx_list)
        contact_mtx = np.concatenate(contact_ensemble, axis=1) if len(contact_ensemble) > 1 else contact_ensemble[0]
        print('contact_mtx shape:', contact_mtx.shape)
        contact_tensor = torch.tensor(contact_mtx, dtype=torch.long).flatten()

        print('contact_tensor shape:', contact_tensor.shape)
        encoded_tensor = F.one_hot(contact_tensor, num_classes=len(chemdata.num2aa))
        print('one-hot encoded_tensor shape:', encoded_tensor.shape)

        encoded_tensor = encoded_tensor.view(contact_mtx.shape[0], -1, encoded_tensor.shape[1])
        print('encoded_tensor reshape:', encoded_tensor.shape)
        encoded_tensor = torch.sum(encoded_tensor, dim=1)
        encoded_tensor = torch.log(encoded_tensor + 1) # apply log(e) to elements
        protein_embedding = torch.cat((embed_mtx, encoded_tensor), dim=1)
        print('protein_embedding shape:', protein_embedding.shape)

        #save to file
        res_idx = [res.id[1] for res in res_list]
        res_name = [res.get_resname() for res in res_list]
        chain_label = [res.get_parent().id for res in res_list]
        data_to_save = {
            'res_idx': res_idx,
            'res_name': res_name,
            'chain_label': chain_label,
            'embedding': protein_embedding
        }
        save_file = os.path.join(save_dir, f'{apo_pdb}.pt')
        torch.save(data_to_save, save_file)