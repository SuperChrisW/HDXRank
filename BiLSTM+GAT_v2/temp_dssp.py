## output dssp (RSA) for each amino acid  ##
from RUN_FIRST import pdb_clean
from output_dssp import output_dssp
import pandas as pd
import os
import warnings

warnings.filterwarnings("ignore")

root_dir = '/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/HDX_MS_dataset/database_collection/feature/AF_Go_complex'
pdb_dir = f'{root_dir}/structure/beta-gamma'
dssp_dir = f'{root_dir}/dssp_files'
if os.path.isdir(dssp_dir) == False:
    os.mkdir(dssp_dir)
fail_list = []


for i in range(1, 101):
    pdb_fname = f'model_{i}' 
    print('start processing', f'{pdb_fname}.pdb')

    if os.path.isfile(f'{root_dir}/dssp_files/{pdb_fname}.dssp.txt'):
        continue
    if os.path.isfile(f'{pdb_dir}/{pdb_fname}.pdb'):
        #try:
        pdb_clean(f'{pdb_dir}/{pdb_fname}.pdb', f'{pdb_dir}/{pdb_fname}_clean.pdb')
        output_dssp(f'{pdb_fname}_C', f'{pdb_dir}/{pdb_fname}_clean.pdb', 'C', f'{root_dir}/dssp_files/beta-gamma')
        #except Exception as e:
        #    fail_list.append(pdb_fname)
        #    print(e)
        #    continue
        break
    else:
        print(f'{pdb_dir}/{pdb_fname}.pdb does not exist')
print(fail_list)