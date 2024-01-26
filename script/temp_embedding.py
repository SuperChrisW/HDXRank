from MSA_embedding import generate_embedding
"""
Created on Wed Jan. 24 2024
@author: Liyao
Assemble all features into one embedding file
"""
import os
import pandas as pd

root_dir = '/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/HDX_MS_dataset/database_collection/feature/CSN3-CSN8'

if __name__ == "__main__":
    dssp_dir = os.path.join(root_dir, 'dssp_files')
    hhm_dir = os.path.join(root_dir, 'hhm_files')
    rigidity_dir = os.path.join(root_dir, 'rigidity_files')
    save_dir = os.path.join(root_dir, 'embedding_files')
    proflex_dir = os.path.join(root_dir, 'proflex_files')
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)

    for file in os.listdir(proflex_dir):
        if 'decomp' in file: continue
        if file.split('.')[0] == '': continue

        parts = file.split('_')
        PDB = parts[:-3]
        PDB = '_'.join(PDB)
        if os.path.isfile(f'{save_dir}/{PDB}.embedding.txt'):
            continue

        try:
            print('processing:', PDB)
            dssp_file = os.path.join(dssp_dir, f'{PDB}.dssp.txt')
            hhm_file = os.path.join(hhm_dir, f'CSN_complex_C.hhm')
            rigidity_file = os.path.join(rigidity_dir, f'rigidity_{PDB}.txt')

            if os.path.isfile(dssp_file) and os.path.isfile(hhm_file) and os.path.isfile(rigidity_file):
                generate_embedding(hhm_file, dssp_file, rigidity_file, f'{save_dir}/{PDB}.embedding.txt', 'A')
            else:
                print(f'file not exist at ',PDB)

        except Exception as e:
            print(f'error {e} at ',PDB)
            continue