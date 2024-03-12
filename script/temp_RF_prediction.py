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
target_HDX_file = '/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/HDX_MS_dataset/database_collection/PXD019367/apo_revised_modified.xlsx'
apo_df = pd.read_excel(target_HDX_file, sheet_name='Sheet1')
apo_df = apo_df[apo_df['protein'] == 'GoA_alpha']
print(apo_df.shape)

root_dir = "/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/HDX_MS_dataset/database_collection/feature/AF_Go_complex"
dssp_dir = os.path.join(root_dir, 'dssp_files/alpha-beta')
hhm_dir = os.path.join(root_dir, 'hhm_files')
rigidity_dir = os.path.join(root_dir, 'rigidity_files/alpha-beta')
proflex_dir = os.path.join(root_dir, 'proflex_files/alpha-beta')

RF_model_path = '/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/HDX_MS_dataset/database_collection/feature/model_files/RF_240123.joblib'
LSTM_model_path = '/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/HDX_MS_dataset/database_collection/feature/model_files/model_weights.pth'
##################################### initial setting ##################################### 

def prepare_data(file, correction_value=0):
    example_embedding = f'{embedding_dir}/{file}'
    if os.path.isfile(example_embedding) == False or os.path.isfile(target_HDX_file) == False:
        print(example_embedding, target_HDX_file)
        return [], []
    input_array, truth_array, start_pos, end_pos, log_t = seq_embedding(target_HDX_file, example_embedding, protein, state, mismatch_report=True, correction_value=correction_value)
    return input_array, truth_array, start_pos, end_pos

def model_prediction(model, embedding_dir, chain_id= 'A', correction_value=0):
    for file in os.listdir(embedding_dir):
        PDB = file.split('.')[0]
        if PDB == '':
            continue
        if not PDB == 'AF_Go_complex_revised_C':
           continue
        parts = PDB.split('_')
        print('processing:', PDB)

        fname = '_'.join(parts[:-1])

        binding_mask = []
        input_array = []
        truth_array = []
        input_array, truth_array, start_pos, end_pos = prepare_data(file, correction_value=correction_value)
        truth_array = np.array(truth_array)/100
        #x = input_array.mean(axis=1)
        x = input_array
        y = truth_array

        print('input shape:', input_array.shape, y.shape)

        pdb_file = f'{root_dir}/structure/alpha-beta/{fname}.pdb'
        chains = read_PDB(PDB, pdb_file)
        Chain_dict = bindingsite_extract(chains, dist_cutoff=4)
        print(Chain_dict)
        for chain_id, sites in Chain_dict.items():
            if chain_id == chain_id:
                min_site = min(sites)
                max_site = max(sites)
        for i, (start, end) in enumerate(zip(start_pos, end_pos)):
            if (start >= min_site and start <= max_site) or (end >= min_site and end <= max_site):
                binding_mask.append(1)
            else:
                binding_mask.append(0)
        binding_mask = np.array(binding_mask)
        x = x[binding_mask == 1, :, :]
        y = y[binding_mask == 1]

        print('filtered by binding site:', x.shape, y.shape)
        if y.shape[0] < 3:
            continue
        model.eval()
        x = torch.tensor(x, dtype=torch.float32)
        x = x.unsqueeze(1)
        with torch.no_grad():
            y_pred = model(x)
        y_pred = y_pred.squeeze().numpy()

        '''
        plt.figure()
        plt.scatter(y[binding_mask == 1], y_pred[binding_mask == 1], c='r', label='binding site')
        plt.scatter(y[binding_mask == 0], y_pred[binding_mask == 0], c='b', label='non-binding site')
        plt.legend()
        plt.xlabel('Experimental HDX')
        plt.ylabel('Predicted HDX')
        plt.title(f'{fname}')
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.show()
        '''

        mse_value = np.mean((y - y_pred) ** 2)
        rmse_value = np.sqrt(mse_value)
        eval_results.append(rmse_value)
        PDB_list.append(fname)

    return eval_results, PDB_list


if __name__ == "__main__":

    predict_df = pd.DataFrame()
    count = 0
    x_array = np.empty((0, 30))
    PDB_list = []
    y_pred_list = []
    eval_results = []
    LSTM_model = CustomModel()
    LSTM_model.load_state_dict(torch.load(LSTM_model_path, map_location=torch.device('cpu')))


    embedding_dir = os.path.join(root_dir, 'embedding_files/alpha')
    state = 'Alone'
    protein = 'GoA_alpha'
    eval_results, PDB_list = model_prediction(LSTM_model, embedding_dir, chain_id='A')
    results_alpha = {pdb: eval for pdb, eval in zip(PDB_list, eval_results)}

    embedding_dir = os.path.join(root_dir, 'embedding_files/beta')
    state = 'Alone'
    protein = 'GoA-beta'
    eval_results, PDB_list = model_prediction(LSTM_model, embedding_dir, chain_id='C', correction_value=-1)
    results_beta = {pdb: eval for pdb, eval in zip(PDB_list, eval_results)}
    #print(results_beta)

    final_results = {}
    for pdb, eval_alpha in results_alpha.items():
        eval_beta = results_beta.get(pdb)   
        final_results[pdb] = eval_alpha * eval_beta

    sorted_final_results_asc = {pdb: score for pdb, score in sorted(final_results.items(), key=lambda item: item[1])}

    for rank, (name, val) in enumerate(sorted_final_results_asc.items()):
        print(f"Rank {rank}: PDB = {name}, RMSE = {val}")


