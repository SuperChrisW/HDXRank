### load LSTM model ###
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
root_dir = "/home/lwang/AI-HDX-main/HDX_MS_dataset/complexpair_dataset"
HDX_summary_file = f'{root_dir}/merged_complex_pair.xlsx'

dssp_dir = os.path.join(root_dir, 'dssp_files')
hhm_dir = os.path.join(root_dir, 'hhm_files')
rigidity_dir = os.path.join(root_dir, 'rigidity_files')
embedding_dir = os.path.join(root_dir, 'embedding_files')
proflex_dir = os.path.join(root_dir, 'proflex_files')
pepGraph_dir = os.path.join(root_dir, 'pepGraph_files')

model_path = '/home/lwang/models/HDX_LSTM/results/240219/model.pt'
save_dir = '/home/lwang/models/HDX_LSTM/results/240219'
##################################### initial setting ##################################### 

if __name__ == "__main__":
    predict_df = pd.DataFrame()
    count = 0
    x_array = np.empty((0, 30))
    PDB_list = []
    y_pred_list = []
    r2_list = []

    device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')
    print(device)
    model = torch.load(model_path, map_location=device)

    apo_df = pd.read_excel(HDX_summary_file, sheet_name='Sheet1')
    apo_df = apo_df.dropna(subset=['chain_identifier'])

    count = 0
    x = np.empty((0, 47, 30))
    y = np.empty(0)
    for index, row in apo_df.iterrows():
        #x = np.empty((0, 47, 30))
        #y = np.empty(0)
        folder = row['database_id'].strip()

        print('processing', folder)
        uni_id = row['match_uni'].strip()
        state = row['state'].strip()
        apo_pdb = row['apo_identifier'].strip().split('.')[0]
        chain_id = row['chain_identifier'].strip()

        target_HDX_file = f'{root_dir}/../{folder}/apo_revised_modified.xlsx'
        target_embedding_file = f'{root_dir}/embedding_files/{apo_pdb}_{chain_id}.embedding.txt'
        if os.path.isfile(target_HDX_file) == False or os.path.isfile(target_embedding_file) == False:
            continue

        input_array, truth_array, start_pos, end_pos, log_t = seq_embedding(target_HDX_file, target_embedding_file, state, mismatch_report=False)
        if input_array.shape[0] == 0:
            continue
        truth_array = np.array(truth_array)/100
        #x = input_array
        #y = truth_array

        x = np.concatenate((x, input_array), axis=0)
        y = np.concatenate((y, truth_array), axis=0)
        print(x.shape, y.shape)

    LSTM_model.eval()
    x = torch.tensor(x, dtype=torch.float32)
    #x = x.permute(0, 2, 1)        
    x = x.unsqueeze(1)
    with torch.no_grad():
        y_pred = LSTM_model(x)
    y_pred = y_pred.squeeze().numpy()
    print(y_pred.shape)

    #### sort prediction by start position ####
    #sorted_lists = sorted(zip(start_pos, end_pos, log_t, y_pred, y), key=lambda x: x[0])
    #start_pos, end_pos, log_t, y_pred, y = zip(*sorted_lists)


    ### plot y_pred vs y ###
    plt.figure()
    plt.scatter(y, y_pred)
    plt.xlabel('True Values [y_test]')
    plt.ylabel('Predictions [y_pred]')
    plt.xlim([0,1])
    plt.ylim([0,1])
    plt.title('True vs Predicted Values')
    plt.show()

    #for index, (start, end, time) in enumerate(zip(start_pos, end_pos, log_t)):
    #    pass
