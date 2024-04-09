### add H to pdb file, then run FIRST ###
import os
import shutil
from RUN_FIRST import RUN_FIRST, copyfiles, pdb_clean, hbplus
import pandas as pd
import os
import warnings


root_dir = '/home/lwang/models/HDX_LSTM/data/COVID_SPIKE'
proflex_dir = f'/home/lwang/AI-HDX-main/ProteinComplex_HDX_prediction/prog/ProFlex-master/proflex'
hbplus_dir = '/home/lwang/AI-HDX-main/ProteinComplex_HDX_prediction/prog/hbplus'

pdb_dir = f'{root_dir}/structure'
save_dir = f'{root_dir}/proflex_files'

if os.path.isdir(save_dir) == False:
    os.mkdir(save_dir)

fname = 'COVID_record.xlsx'
apo_df = pd.read_excel(f'{root_dir}/{fname}', sheet_name='Sheet1')
apo_df = apo_df.dropna(subset=['chain_identifier'])

fail_list = []
for index, row in apo_df.iterrows():
    pdb_fname = row['apo_identifier'].strip().split('.')[0]
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