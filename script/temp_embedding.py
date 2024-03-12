from MSA_embedding import generate_embedding
"""
Created on Wed Jan. 24 2024
@author: Liyao
Assemble all features into one embedding file
"""
import os
import pandas as pd

root_dir = '/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/HDX_MS_dataset/database_collection/feature/AF_Go_complex'

if __name__ == "__main__":
    dssp_dir = os.path.join(root_dir, 'dssp_files/alpha-beta')
    hhm_dir = os.path.join(root_dir, 'hhm_files')
    rigidity_dir = os.path.join(root_dir, 'rigidity_files/alpha-beta')
    save_dir = os.path.join(root_dir, 'embedding_files/alpha-beta')
    proflex_dir = os.path.join(root_dir, 'proflex_files/alpha-beta')
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)
    
    chain_id = 'C'
    for file in os.listdir(proflex_dir):
        if 'decomp' in file: continue
        if file.split('.')[0] == '': continue

        parts = file.split('_')
        PDB = parts[:-3]
        PDB = '_'.join(PDB)
        if os.path.isfile(f'{save_dir}/{PDB}_{chain_id}.embedding.txt'):
            continue

#try:
        print('processing:', PDB)
        dssp_file = os.path.join(dssp_dir, f'{PDB}_{chain_id}.dssp.txt')
        hhm_file = os.path.join(hhm_dir, f'AF_Go_complex_revised_{chain_id}.hhm')
        rigidity_file = os.path.join(rigidity_dir, f'rigidity_{PDB}.txt')

        if os.path.isfile(dssp_file) and os.path.isfile(hhm_file) and os.path.isfile(rigidity_file):
            generate_embedding(hhm_file, dssp_file, rigidity_file, f'{save_dir}/{PDB}_{chain_id}.embedding.txt', chain_id)
        else:
            print(f'file not exist at ',PDB)
    '''
        except Exception as e:
            print(f'error {e} at ',PDB)
            continue
        '''