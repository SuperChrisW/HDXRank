### finetune the LSTM model toward the target protein ###
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
target_embedding_file = '/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/HDX_MS_dataset/database_collection/feature/CSN3-CSN8/embedding_files/6CSN8.fasta'
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

LSTM_model_path = '/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/HDX_MS_dataset/database_collection/feature/model_files/model_weights.pth'
##################################### initial setting ##################################### 

def prepare_data(target_HDX_file, target_embedding_file):
    if os.path.isfile(target_embedding_file) == False or os.path.isfile(target_HDX_file) == False:
        print(target_embedding_file, target_HDX_file)
        return [], []
    input_array, truth_array, start_pos, end_pos, log_t = seq_embedding(target_HDX_file, target_embedding_file, state, mismatch_report=False)
    return input_array, truth_array, start_pos, end_pos

if __name__ == "__main__":
    LSTM_model = CustomModel()
    LSTM_model.load_state_dict(torch.load(LSTM_model_path, map_location=torch.device('cpu')))