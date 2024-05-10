"""
Created on Wed Jan. 24 2024

@author: liyao
convert proflex output to ASSURR flexibility index
Usage: python temp_flexibility.py
"""

import os
from prot_rigidity import ASSURR_flexibility

proflex_dir = '/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/HDX_MS_dataset/database_collection/feature/AF_Go_complex/proflex_files/alpha-beta'
dssp_dir = '/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/HDX_MS_dataset/database_collection/feature/AF_Go_complex/dssp_files/alpha-beta'
save_dir = '/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/HDX_MS_dataset/database_collection/feature/AF_Go_complex/rigidity_files/alpha-beta'


if __name__ == '__main__':
    if not os.path.exists(proflex_dir) and not os.path.exists(dssp_dir):
        print('cannot find the folder')
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    for file in os.listdir(proflex_dir):
        if 'decomp' in file: continue
        if file.split('.')[0] == '': continue

        parts = file.split('_')
        PDB = parts[:-3]
        PDB = '_'.join(PDB)
        if os.path.exists(f'{save_dir}/rigidity_{PDB}.txt'): continue

        #skip_list = ['6h9v', 'CSN_complex', 'CD47_BRIL', 'DAT+DA_complex', '1t3d_hexamer']

        dssp_file = []
        dssp_fname = [f'{PDB}_A.dssp.txt', f'{PDB}_C.dssp.txt']
        for fname in dssp_fname:
            dssp_file.append(f'{dssp_dir}/{fname}')

        ASSURR_flexibility(PDB, dssp_file, proflex_dir, save_dir)