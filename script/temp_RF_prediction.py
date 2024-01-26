"""
Created on Wed Jan. 24, 2024
@author: Liyao
Load the RF model and make prediction for target protein
"""
import pandas as pd
import os
from polyR_model import seq_embedding
import numpy as np
import joblib
from sklearn.metrics import r2_score, mean_squared_error
from scipy.stats import pearsonr, zscore
import matplotlib.pyplot as plt
from binding_site_extraction import read_PDB, bindingsite_extract

import torch
from torch_lstm import CustomModel
##################################### initial setting ##################################### 
target_HDX_file = '/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/HDX_MS_dataset/database_collection/PXD013001/apo_revised_modified.xlsx'
apo_df = pd.read_excel(target_HDX_file, sheet_name='Sheet1')
apo_df = apo_df[apo_df['state'] == 'CSN3']
print(apo_df.shape)
state = 'CSN3'

root_dir = "/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/HDX_MS_dataset/database_collection/feature/CSN3-CSN8"
dssp_dir = os.path.join(root_dir, 'dssp_files')
hhm_dir = os.path.join(root_dir, 'hhm_files')
rigidity_dir = os.path.join(root_dir, 'rigidity_files')
embedding_dir = os.path.join(root_dir, 'embedding_files')
proflex_dir = os.path.join(root_dir, 'proflex_files')

RF_model_path = '/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/HDX_MS_dataset/database_collection/feature/model_files/RF_240123.joblib'
LSTM_model_path = '/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/HDX_MS_dataset/database_collection/feature/model_files/model_240125.pt'
##################################### initial setting ##################################### 

def prepare_data():
    example_embedding = f'{embedding_dir}/{file}'
    if os.path.isfile(example_embedding) == False or os.path.isfile(target_HDX_file) == False:
        print(example_embedding, target_HDX_file)
        return [], []
    input_array, truth_array, start_pos, end_pos = seq_embedding(target_HDX_file, example_embedding, state, mismatch_report=False)
    return input_array, truth_array, start_pos, end_pos


if __name__ == "__main__":
    #predict_df = apo_df[['state', 'sequence', 'exposure', 'start', 'end', '%d']]
    predict_df = pd.DataFrame()
    count = 0
    x_array = np.empty((0, 30))
    PDB_list = []
    y_pred_list = []
    r2_list = []
    #random_forest_model = joblib.load(model_path)
    LSTM_model = torch.load(LSTM_model_path, map_location=torch.device('cpu'))

    for file in os.listdir(embedding_dir):
        PDB = file.split('.')[0]
        print('processing:', PDB)
        pdb_file = f'{root_dir}/structure/{PDB}.pdb'
        chains = read_PDB(PDB, pdb_file)
        Chain_dict = bindingsite_extract(chains, dist_cutoff=3.65)

        for chain_id, sites in Chain_dict.items():
            if chain_id == 'A':
                min_site = min(sites)
                max_site = max(sites)
        
        binding_mask = []
        input_array = []
        truth_array = []
        input_array, truth_array, start_pos, end_pos = prepare_data()
        truth_array = np.array(truth_array)/100
        #x = input_array.mean(axis=1)
        x = input_array
        y = truth_array
        #x = x + 1e-10
        print(input_array.shape, y.shape)

        for i, (start, end) in enumerate(zip(start_pos, end_pos)):
            if (start >= min_site and start <= max_site) or (end >= min_site and end <= max_site):
                binding_mask.append(1)
            else:
                binding_mask.append(0)
        binding_mask = np.array(binding_mask)
        x = x[binding_mask == 1, :, :]
        y = y[binding_mask == 1]

        LSTM_model.eval()
        x = torch.tensor(x, dtype=torch.float32)
        x = x.permute(0, 2, 1)        
        x = x.unsqueeze(1)
        with torch.no_grad():
            y_pred = LSTM_model(x)
        y_pred = y_pred.squeeze().numpy()

        y = zscore(y)
        y_pred = zscore(y_pred)
        r2 = r2_score(y, y_pred)
        r2_list.append(r2)
        PDB_list.append(PDB)
    
    prediction_list = [(r2, PDB) for r2, PDB in sorted(zip(r2_list, PDB_list), reverse=True)]

    for rank, (r2, PDB) in enumerate(prediction_list, start=1):
        print(f"Rank {rank}: PDB = {PDB}, r2 Score = {r2}")


'''
    y_pred = random_forest_model.predict(x_array)
    y_pred = np.array(y_pred)
    y_pred = y_pred.reshape(input_array.shape[0], -1)
    predict_df = pd.concat([predict_df, pd.DataFrame(y_pred)], axis=1)
    predict_df.columns = PDB_list
    predict_df['%d'] = y
    predict_df.to_excel(f'{root_dir}/CSN3_prediction.xlsx', index=False)
'''




