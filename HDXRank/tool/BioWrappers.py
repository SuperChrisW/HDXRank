# this file is from DeepRank

import os
import numpy as np
import pandas as pd
from Bio import PDB
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.HSExposure import HSExposureCA

import warnings
from Bio import BiopythonWarning

import biotite.structure as struc
import biotite.structure.io as strucio
import numpy as np

with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonWarning)

from tool.predataprocess import get_pdb_seq
from tool.RUN_FIRST import RUN_FIRST, copyfiles, pdb_clean, hbplus
from tool.prot_rigidity import ASSURR_flexibility

def get_bio_model(pdbfile):
    """Get the model

    Args:
        pdbfile (str): pdbfile

    Returns:
        [type]: Bio object
    """
    parser = PDBParser(QUIET = True)
    structure = parser.get_structure('_tmp', pdbfile)
    return structure[0]

def get_hse(model, chain='A'):
    """Get the hydrogen surface exposure

    Args:
        model (bio model): model of the strucrture

    Returns:
        dict: hse data
    """

    hse = HSExposureCA(model)
    data = {}
    hse_mtx = []
    index_range = []
    for k in list(hse.keys()):
        if not k[0] == chain:
            continue
        new_key = (k[0], k[1][1])
        index_range.append(k[1][1])

        x = hse[k]
        if x[2] is None:
            x = list(x)
            x[2] = 0.0
            x = tuple(x)

        data[new_key] = x
        hse_mtx.append(list(x))
    return data

def biotite_SASA(pdb_file, chain_id, HDXparser = None):
    '''
    biotite.structure.sasa(array, probe_radius=1.4, atom_filter=None, ignore_ions=True, point_number=1000, point_distr='Fibonacci', vdw_radii='ProtOr')
    '''
    model = struc.io.load_structure(pdb_file)
    chain_mask = (model.chain_id == chain_id)

    atom_sasa = struc.sasa(model, atom_filter = chain_mask, vdw_radii="Single")
    res_sasa = struc.apply_residue_wise(model, atom_sasa, np.sum)
    return res_sasa[~np.isnan(res_sasa)]

'''
def output_fasta():
    ### output fasta file for hhblits ###

    root_dir = '/home/lwang/models/HDX_LSTM/data/test_set'
    fname = 'merged_data.xlsx'
    apo_df = pd.read_excel(f'{root_dir}/{fname}', sheet_name='test_set')
    apo_df = apo_df.dropna(subset=['chain_identifier'])
    fail_list = []

    for index, row in apo_df.iterrows():
        pdb_fname = row['apo_identifier'].strip().split('.')[0]
        #chain_id = row['chain_identifier'].strip()
        #uni_id = row['match_uni'].strip()
        print('start processing', f'{pdb_fname}')        
        if os.path.isfile(f'{root_dir}/structure/{pdb_fname}.pdb'):
            try:
                parser = PDB.PDBParser()
                filepath = f'{root_dir}/structure/{pdb_fname}.pdb'
                structure = parser.get_structure('AF_structure', filepath)
                chain_ids = list(structure.get_chains())
                for chain_id in chain_ids:
                    chain_id = chain_id.get_id()
                    sequence, residue_indices = get_pdb_seq(structure, chain_id)
                    if os.path.isdir(f'{root_dir}/fasta_files') == False:
                        os.mkdir(f'{root_dir}/fasta_files')
                    with open(f'{root_dir}/fasta_files/{pdb_fname}_{chain_id}.fasta', 'w') as f:
                        f.write(f'>{pdb_fname}_{chain_id}\n')
                        f.write(sequence)
            except Exception as e:
                print(e)
                fail_list.append(pdb_fname+'_'+chain_id)
                continue
        else:
            print(f'{root_dir}/structure/{pdb_fname}.pdb does not exist')
    print(fail_list)
'''

def output_proflex():
    ### add H to pdb file, then run FIRST ###
    ### usually need to run separately in terminal ###

    root_dir = '/home/lwang/models/HDX_LSTM/data/test_set'
    proflex_dir = f'/home/lwang/AI-HDX-main/ProteinComplex_HDX_prediction/prog/ProFlex-master/proflex'
    hbplus_dir = '/home/lwang/AI-HDX-main/ProteinComplex_HDX_prediction/prog/hbplus'

    pdb_dir = f'{root_dir}/structure'
    save_dir = f'{root_dir}/proflex_files'

    if os.path.isdir(save_dir) == False:
        os.mkdir(save_dir)

    fname = 'merged_data.xlsx'
    apo_df = pd.read_excel(f'{root_dir}/{fname}', sheet_name='test_set')
    apo_df = apo_df.dropna(subset=['chain_identifier'])

    for index, row in apo_df.iterrows():
        pdb_fname = row['apo_identifier'].strip().split('.')[0]
        #chain_id = row['chain_identifier'].strip()
        #uni_id = row['match_uni'].strip()

        pdb_fpath = f'{root_dir}/structure/{pdb_fname}.pdb'

        proflex_name = f'{pdb_fname}_clean_Hplus_proflexdataset'
        if os.path.isfile(f'{save_dir}/{proflex_name}'):
            continue           
        
        if os.path.isfile(pdb_fpath):
            try:
                print('start processing', f'{pdb_fname}') 
                pdb_folder = f'{root_dir}/structure/'
                if pdb_clean(pdb_fpath, f'{hbplus_dir}/{pdb_fname}_clean.pdb') == False:
                    raise Exception('pdb_clean failed')
                hbplus(hbplus_dir, f'{pdb_fname}_clean', pdb_folder)

                pdb_fname = pdb_fname+'_clean_Hplus'       
                pdb_clean(f'{pdb_folder}/{pdb_fname}.pdb', f'{proflex_dir}/{pdb_fname}.pdb')
                RUN_FIRST(proflex_dir, f'{pdb_fname}.pdb', '-h')
                copyfiles(pdb_fname, proflex_dir, save_dir)

            except Exception as e:
                fail_list.append(pdb_fname)
                print(e)
                continue
    print(fail_list)

def output_rigidity():
    root_dir = '/home/lwang/models/HDX_LSTM/data/test_set'
    proflex_dir = f'{root_dir}/proflex_files'
    dssp_dir = f'{root_dir}/dssp_files'
    save_dir = f'{root_dir}/rigidity_files'
    if os.path.isdir(save_dir) == False:
        os.mkdir(save_dir)

    fname = 'merged_data.xlsx'
    apo_df = pd.read_excel(f'{root_dir}/{fname}', sheet_name='test_set')
    apo_df = apo_df.dropna(subset=['chain_identifier'])
    fail_list = []

    for index, row in apo_df.iterrows():
        pdb_fname = row['apo_identifier'].strip().split('.')[0]
        print('start processing', f'{pdb_fname}')        
        if os.path.isfile(f'{root_dir}/structure/{pdb_fname}.pdb'):
            parser = PDB.PDBParser()
            filepath = f'{root_dir}/structure/{pdb_fname}.pdb'
            structure = parser.get_structure('AF_structure', filepath)
            chain_ids = list([chain.id for chain in structure.get_chains()])
            dssp_file = [f'{dssp_dir}/{pdb_fname}_{chain_id}.dssp.txt' for chain_id in chain_ids]
            ASSURR_flexibility(pdb_fname, dssp_file, proflex_dir, save_dir)

