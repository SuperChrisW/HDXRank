from RUN_FIRST import pdb_clean
from output_dssp import output_dssp
import pandas as pd
import os
import warnings

warnings.filterwarnings("ignore")

root_dir = '/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/HDX_MS_dataset/database_collection/feature'
pdb_dir = f'{root_dir}/structure'
dssp_dir = f'{root_dir}/dssp_files'
if os.path.isdir(dssp_dir) == False:
    os.mkdir(dssp_dir)

fname = 'merged_apo.xlsx'
apo_df = pd.read_excel(f'{root_dir}/../{fname}', sheet_name='Sheet1')
apo_df = apo_df.dropna(subset=['chain_identifier'])
fail_list = []

for index, row in apo_df.iterrows():
    pdb_fname = row['apo_identifier'].strip().split('.')[0]
    chain_id = row['chain_identifier'].strip()
    uni_id = row['match_uni'].strip()
    
    print('start processing', f'{pdb_fname}_{chain_id}')        
    if os.path.isfile(f'{root_dir}/dssp_files/{pdb_fname}_{chain_id}.dssp.txt'):
        continue
    if os.path.isfile(f'{root_dir}/structure/{uni_id}/{pdb_fname}.pdb'):
        try:
            pdb_clean(f'{root_dir}/structure/{uni_id}/{pdb_fname}.pdb', f'{root_dir}/structure/{uni_id}/{pdb_fname}_clean.pdb')
            print(f'{root_dir}/structure/{uni_id}/{pdb_fname}_clean.pdb')
            output_dssp(f'{pdb_fname}_{chain_id}', f'{root_dir}/structure/{uni_id}/{pdb_fname}_clean.pdb', chain_id, f'{root_dir}/dssp_files')
        except Exception as e:
            fail_list.append(pdb_fname)
            print(e)
            continue
    else:
        print(f'{root_dir}/structure/{uni_id}/{pdb_fname}.pdb does not exist')
print(fail_list)