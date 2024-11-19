# pre-setting and load prediction
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.stats import spearmanr
import pickle

def get_true_diff(HDX_fpath, apo_states, complex_states, cluster_id, timepoints):
    def get_weighted_uptake(HDX_df, protein, state, correction):
        temp_HDX_df = HDX_df[(HDX_df['state']==state) & (HDX_df['protein']==protein)]
        if cluster_id == 0:
            temp_HDX_df = temp_HDX_df[temp_HDX_df['log_t']<timepoints[0]]
        elif cluster_id == 1:
            temp_HDX_df = temp_HDX_df[(temp_HDX_df['log_t']<timepoints[1]) & (temp_HDX_df['log_t']>=timepoints[0])]
        elif cluster_id == 2:
            temp_HDX_df = temp_HDX_df[(temp_HDX_df['log_t']>=timepoints[1])]
        else:
            raise ValueError('Invalid cluster_id')
        
        temp_HDX_df = temp_HDX_df.sort_values(by=['start', 'end'], ascending=[True, True])

        exposures = temp_HDX_df['exposure'].unique()
        unweighted_RFU = {time:{} for time in exposures}
        for time in exposures:
            for index, row in temp_HDX_df[temp_HDX_df['exposure']==time].iterrows():
                #chain_id = int(row['Chain'])
                unweighted_RFU[time][f'{row["start"]+correction}-{row["end"]+correction}'] = row['RFU']

        grouped = temp_HDX_df.groupby(['start', 'end'])
        weighted_uptake = grouped['RFU'].mean().to_numpy()/100
        unique_pairs = grouped.groups.keys()
        x_label = [f'{start+correction}-{end+correction}' for start, end in unique_pairs]

        return weighted_uptake, x_label, unweighted_RFU
    
    HDX_df = pd.read_excel(HDX_fpath)
    true_apo, true_complex = {}, {}
    apo_mtx, complex_mtx = {}, {}

    protein, state, correction = apo_states[0], apo_states[1], apo_states[2]
    uptake, label, U_apo_mtx = get_weighted_uptake(HDX_df, protein, state, correction)
    for l, u in zip(label, uptake):
        true_apo[l] = u
    apo_mtx.update(U_apo_mtx)

    protein, state, correction = complex_states[0], complex_states[1], complex_states[2]
    uptake, label, U_complex_mtx = get_weighted_uptake(HDX_df, protein, state, correction)
    for l, u in zip(label, uptake):
        true_complex[l] = u
    complex_mtx.update(U_complex_mtx)

    true_diff = {}
    diff_mtx = {}
    for key in true_apo.keys():
        if key in true_complex:
            true_diff[key] = true_complex[key] - true_apo[key]
            diff_mtx[key] = {t: complex_mtx[t][key] - apo_mtx[t][key] for t in apo_mtx if key in apo_mtx[t] and key in complex_mtx[t]}
    
    print('Common peptides num:', len(true_diff.keys()))
    return true_diff, diff_mtx

def plot_HDX_diff(true_diff, diff_mtx, size=(10, 6)):
    x_labels = list(true_diff.keys())
    diff = np.array(list(true_diff.values()))
    diff_neg = diff[diff<0]
    mean_diff = np.mean(diff_neg)
    hdx_epitope_id = np.where(diff<mean_diff)[0]
    hdx_epitope_pep = [x_labels[i] for i in hdx_epitope_id]

    all_times = set()
    for diffs in diff_mtx.values():
        all_times.update(diffs.keys())
    sorted_times = sorted(all_times, key=lambda x: float(x))

    return hdx_epitope_id, hdx_epitope_pep, sorted_times

def root_mean_square_error(y_true, y_pred, error_limit=1):
    return np.mean((((y_true - y_pred)/error_limit) ** 2))

def prepare_data(complex_batch, apo_states, complex_states, hdx_true_diffs, hdx_epitope_peps=None):
    truth = []
    pred = []
    if complex_batch not in average_df['Batch'].values:
        return None, None
    if hdx_epitope_peps is None:
        hdx_epitope_peps = [list(hdx_dict.keys()) for hdx_dict in hdx_true_diffs]
    for apo_state, complex_state, epitope_peps, hdx_dict in zip(apo_states, complex_states, hdx_epitope_peps, hdx_true_diffs):
        apo_batch = apo_state[-1]
        #complex_batch = complex_state[-1]
        apo_df = average_df[average_df['Batch']==apo_batch]
        complex_df = average_df[average_df['Batch']==complex_batch]

        for pep in epitope_peps:
            if pep in apo_df['Range'].values and pep in complex_df['Range'].values:
                pred_diff = complex_df.loc[complex_df['Range']==pep, 'Y_Pred'].values[0] - apo_df.loc[apo_df['Range']==pep, 'Y_Pred'].values[0]
                true_diff = hdx_dict[pep]
                truth.append(true_diff)
                pred.append(pred_diff)
    return np.array(truth), np.array(pred)

def parse_dockq_results(pkl_file):
    with open(pkl_file, 'rb') as f:
        dockq_results = pickle.load(f)
    data = []
    model_list = [f'MODEL_{i}_REVISED' for i in range(1, len(dockq_results)+1)]
    for result in dockq_results:
        entry = {
            'lrms': result.get('Lrms', np.nan),
            'irms': result.get('irms', np.nan),
            'fnat': result.get('fnat', np.nan),
            'dockq': result.get('DockQ', np.nan)
        }
        data.append(entry)
    df = pd.DataFrame(data, columns=['lrms', 'irms', 'fnat', 'dockq'])
    df['Batch'] = model_list
    return df

def get_HDX_score(apo_states, complex_states, HDX_fpath, cluster_id, timepoints, N_decoys=1000):
    hdx_epitope_peps = []
    hdx_true_diffs = []
    for apo, complex in zip(apo_states, complex_states):
        true_diff, diff_mtx = get_true_diff(HDX_fpath, apo, complex, cluster_id, timepoints)
        epitope_id, epitope_pep, hdx_times = plot_HDX_diff(true_diff, diff_mtx, size=(12,4))
        hdx_true_diffs.append(true_diff)
        hdx_epitope_peps.append(epitope_pep)

    HDX_scores = []
    complex_batch_list = [f'MODEL_{i}_REVISED' for i in range(1, N_decoys+1)]

    for complex_batch in tqdm(complex_batch_list):
        y_true, y_pred = prepare_data(complex_batch, apo_states, complex_states, hdx_true_diffs, hdx_epitope_peps=hdx_epitope_peps)
        if y_true is None:
            continue
        hdx_score = root_mean_square_error(y_true, y_pred)
        HDX_scores.append(hdx_score)

    hdx_times = [str(t) for t in hdx_times]
    HDX_score_df = pd.DataFrame({'Batch': complex_batch_list, 'HDX_score': HDX_scores, 'HDX_time': ','.join(hdx_times)})
    dockq_df = parse_dockq_results(f'{root_dir}/{protein_name}/eval/{protein_name}_dockq.pkl')
    total_df = pd.merge(HDX_score_df, dockq_df, on='Batch')

    spr_corr, _ = spearmanr(total_df['lrms'], total_df['HDX_score'])
    print(f'Spearman correlation: {spr_corr:.4f}')


root_dir = '/home/lwang/models/HDX_LSTM/data/hdock/structure'
protein_name = '8F7A_ori'
cluster_id = 1
N_decoys = 1000
cluster = f'cluster{cluster_id}_8A_manual_rescale'
HDX_fpath = '/home/lwang/models/HDX_LSTM/data/hdock/HDX_files/PXD037571_revised.xlsx'

dfs = []
epochs = [60, 70, 80]
for i,epoch in enumerate(epochs):
    #fpath = f'/home/lwang/models/HDX_LSTM/data/hdock/prediction/1016/GN56_epoch{epoch}/HDX_pred_GearNet56_{protein_name}_{cluster}_v{i+1}.csv'
    fpath = f'/home/lwang/models/HDX_LSTM/data/hdock/prediction/GN56_epoch{epoch}/HDX_pred_GearNet56_{protein_name}_{cluster}_v0.csv'
    df = pd.read_csv(fpath)
    print(df.shape)
    dfs.append(df)
merged_df = pd.concat(dfs, axis=0)
average_df = merged_df.groupby(['Batch', 'Chain', 'Range', 'Y_True'])['Y_Pred'].mean().reset_index()

summary_df = pd.read_excel('/home/lwang/models/HDX_LSTM/data/hdock/test_data_AF.xlsx', sheet_name='hdock')
summary_df = summary_df.dropna(subset=['match_uni'])
summary_df = summary_df[(summary_df['note']==protein_name[:4])]
grouped = summary_df.groupby(['match_uni'])

apo_states, complex_states = [], []
for name, group in grouped:
    apo_group = group[group['complex_state'] == 'single']
    complex_group = group[group['complex_state'] == 'protein complex']
    apo_states.append((apo_group['protein'].values[0], apo_group['state'].values[0], apo_group['correction_value'].values[0], apo_group['structure_file'].values[0]))
    complex_states.append((complex_group['protein'].values[0], complex_group['state'].values[0], complex_group['correction_value'].values[0], complex_group['structure_file'].values[0]))

for i in range(1,4):
    timepoints = [i, i+1]
    print(f'protein: {protein_name}, Timepoints {timepoints}')
    get_HDX_score(apo_states, complex_states, HDX_fpath, cluster_id, timepoints, N_decoys=N_decoys)




