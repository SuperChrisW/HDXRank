"""
Created on Wed Jan. 24 2024

@author: liyao
convert proflex output to ASSURR flexibility index
Usage: python temp_flexibility.py
"""

import os
from prot_rigidity import ASSURR_flexibility

proflex_dir = '/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/HDX_MS_dataset/database_collection/feature/CSN3-CSN8/proflex_files'
dssp_dir = '/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/HDX_MS_dataset/database_collection/feature/CSN3-CSN8/dssp_files'
save_dir = '/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/HDX_MS_dataset/database_collection/feature/CSN3-CSN8/rigidity_files'


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

        dssp_fname = f'{PDB}.dssp.txt'
        dssp_file = f'{dssp_dir}/{dssp_fname}'

        ASSURR_flexibility(PDB, dssp_file, proflex_dir, save_dir)