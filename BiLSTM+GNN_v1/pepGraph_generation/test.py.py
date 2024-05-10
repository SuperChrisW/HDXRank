import os
import pandas as pd
import torch
from pepGraph_basemodel_training import base_data, test_basemodel
import tensorflow as tf
import random

##################################### initial setting ##################################### 
root_dir = "/home/lwang/models/HDX_LSTM/data/COVID_SPIKE"
HDX_summary_file = f'{root_dir}/COVID_record.xlsx'
hdx_df = pd.read_excel(HDX_summary_file, sheet_name='Sheet1')
hdx_df = hdx_df.dropna(subset=['chain_identifier'])

dssp_dir = os.path.join(root_dir, 'dssp_files')
hhm_dir = os.path.join(root_dir, 'hhm_files')
rigidity_dir = os.path.join(root_dir, 'rigidity_files')
proflex_dir = os.path.join(root_dir, 'proflex_files')
embedding_dir = os.path.join(root_dir, 'embedding_files')
HDX_dir = os.path.join(root_dir, 'HDX_files')

result_dir = '/home/lwang/models/HDX_LSTM/data/COVID_SPIKE'
if not os.path.exists(result_dir):
    os.mkdir(result_dir)
##################################### initial setting ##################################### 
    
device = torch.device('cpu')
print(device)

### data preparation ###
input_apo, truth_apo, input_complex, truth_complex, range_list = base_data(hdx_df, root_dir)

print('length of input_data:', len(input_apo)+len(input_complex))
input_apo = tf.convert_to_tensor(input_apo, dtype=tf.float32)
truth_apo = tf.convert_to_tensor(truth_apo, dtype=tf.float32)
input_complex = tf.convert_to_tensor(input_complex, dtype=tf.float32)
truth_complex = tf.convert_to_tensor(truth_complex, dtype=tf.float32)