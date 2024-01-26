### add H to pdb file, then run FIRST ###
import os
import shutil
from RUN_FIRST import RUN_FIRST, copyfiles, pdb_clean, hbplus
import pandas as pd
import os
import warnings

root_dir = '/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/HDX_MS_dataset/database_collection/feature/CSN3-CSN8'
proflex_dir = f'/Users/liyao/Desktop/Tsuda_Lab/Source_code/ProFlex-master/proflex'
hbplus_dir = '/Users/liyao/Desktop/Tsuda_Lab/Source_code/hbplus'

pdb_dir = f'{root_dir}/structure'
save_dir = f'{root_dir}/proflex_files'
fail_list = []
if os.path.isdir(save_dir) == False:
    os.mkdir(save_dir)
'''
fname = 'merged_apo.xlsx'
apo_df = pd.read_excel(f'{root_dir}/../{fname}', sheet_name='Sheet1')
apo_df = apo_df.dropna(subset=['chain_identifier'])
'''

for i in range(1, 101):
    pdb_fname = f'model_{i}_revised' 
    print('start processing', f'{pdb_fname}.pdb')  
    
    proflex_name = f'{pdb_fname}_clean_Hplus_proflexdataset'
    if os.path.isfile(f'{save_dir}/{proflex_name}'):
        continue           

    pdb_fpath = f'{pdb_dir}/{pdb_fname}.pdb'
    
    if os.path.isfile(pdb_fpath):
        try:
            if pdb_clean(pdb_fpath, f'{hbplus_dir}/{pdb_fname}_clean.pdb') == False:
                raise Exception('pdb_clean failed')
            hbplus(hbplus_dir, f'{pdb_fname}_clean.pdb', pdb_dir)

            pdb_fname = pdb_fname+'_clean_Hplus'       
            pdb_clean(f'{pdb_dir}/{pdb_fname}.pdb', f'{proflex_dir}/{pdb_fname}.pdb')
            RUN_FIRST(proflex_dir, f'{pdb_fname}.pdb', '-h')
            copyfiles(pdb_fname, proflex_dir, save_dir)

        except Exception as e:
            fail_list.append(pdb_fname)
            print(e)
            continue

print(fail_list)