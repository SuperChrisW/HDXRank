from pyhdx import read_dynamx, HDXMeasurement
from pyhdx.plot import peptide_coverage
from pyhdx.process import correct_d_uptake, filter_peptides, apply_control
from pyhdx.models import Coverage
from pyhdx.config import cfg
import proplot as pplt
import pandas as pd
import numpy as np
import os 

# get residue-wise HDX uptake 
from pyhdx.fileIO import csv_to_hdxm, csv_to_dataframe
from pyhdx.fitting import fit_d_uptake
import pickle

root_dir = f'/home/lwang/models/HDX_LSTM/data/Latest_set'
hdx_dir = f'{root_dir}/HDX_files/'
pdb_dir = f'{root_dir}/structure/'
hdxm_dir = f'{root_dir}/HDX_files/hdxm_files'
save_dir = f'{root_dir}/HDX_files/res_D_fit'

save_fail = []
for file in os.listdir(hdxm_dir):
    file_id = int(file.split('.')[0])    
    if os.path.exists(f'{save_dir}/{file_id}.pkl'):
        continue
    
    try:
        df = csv_to_dataframe(f'{hdxm_dir}/{file}')
        sequence = df.attrs["metadata"]['sequence'].replace('x', '')
        hdxm = csv_to_hdxm(f'{hdxm_dir}/{file}')
        
        fit_1 = fit_d_uptake(hdxm, r1=1.0, repeats=20)
        result = (sequence, fit_1.result)
        #assert fit_1.result.shape[-1] == len(sequence)-1

        with open(f'{save_dir}/{file_id}.pkl', 'wb') as f:
            pickle.dump(result, f)
    except Exception as e:
        print(f'Failed in file{file_id}:', e)
        save_fail.append(file)
        continue
print('Failed list:', save_fail)