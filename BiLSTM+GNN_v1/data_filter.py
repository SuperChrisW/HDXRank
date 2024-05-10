from os import remove
import pandas as pd
from data_approximations import linreg
from protein_merge import protein_merge
from simplemodel import vector, drate, removeOverlap, training
import os
import numpy as np

root_dir = "/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/dataset"
manager_file = f"{root_dir}/filtered_results/manager_finished.txt"
manager = open(manager_file, "r")
hdx_total = []
input_total = []
for line in manager:
    line = line.split(':')
    if line[1] != '':
        csv_file = line[0]
        identifier = line[1].strip()
        prediction_time = line[2].strip()
        Uni_ID = line[3].strip()
        correction_value = line[4].strip()
        if correction_value == 'nan':
            correction_value = 0
            pass
            #continue
    else:
        continue
    
    if os.path.exists(f'{root_dir}/{csv_file}.csv'):
        df = pd.read_csv(f'{root_dir}/{csv_file}.csv')
    else:
        continue

    trim_df = df[df['Exposure'].astype(int) == int(prediction_time)]
    trim_df = trim_df[trim_df['State'] == identifier]
    trim_df = trim_df[trim_df['Start']>0]
    trim_df['Start'] = trim_df['Start'].astype(int)+int(correction_value)
    trim_df['End'] = trim_df['End'].astype(int)+int(correction_value)
    
    #maxuptake
    if "MaxUptake" not in trim_df.columns:
        if 'Sequence' in trim_df.columns:
            trim_df['MaxUptake'] = trim_df['Sequence'].apply(lambda x: len(x) - 1 - x.count('P'))
            pass

    trim_df['MaxUptake'] = trim_df['MaxUptake'].astype(int)
    Uptake_columns = ['Uptake', 'Uptake (Da)', '#D']
    back_exchange_rate = 0.7

    for col in Uptake_columns:
        if col in trim_df.columns:
            trim_df['%D'] = trim_df[col] / (trim_df['MaxUptake'] * back_exchange_rate)
            over_range_mask = trim_df['%D'] > 1
            below_range_mask = trim_df['%D'] < 0
            trim_df.loc[over_range_mask, '%D'] = 1
            trim_df.loc[below_range_mask, '%D'] = 0

    trim_df.to_csv(f'{root_dir}/filtered_results/{csv_file}_filtered.csv', index=False)
    print(csv_file, identifier, prediction_time, Uni_ID, correction_value)

#print('training')
#print(len(input_total), len(hdx_total))
#training(input_total, hdx_total)
    
    
 
