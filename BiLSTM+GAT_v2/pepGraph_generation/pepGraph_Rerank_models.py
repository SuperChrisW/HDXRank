# re-rank of Hdock models based on predicted HDX value
# 2024.04.20 by WANG LIYAO

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LinearRegression
from Bio.PDB import PDBParser, Superimposer
from pdb2sql import interface

from pepGraph_utlis import read_PDB, bindingsite_extract

import itertools
import warnings 
warnings.filterwarnings("ignore")

#get the true HDX difference from source file
def get_true_diff(HDX_fpath, apo_states, complex_states):
    def get_weighted_uptake(HDX_df, protein, state, chain):
        temp_HDX_df = HDX_df[(HDX_df['state']==state) & (HDX_df['protein']==protein)]
        temp_HDX_df = temp_HDX_df.sort_values(by=['start', 'end'], ascending=[True, True])
        
        exposures = temp_HDX_df['log_t'].unique()
        unweighted_RFU = []
        for time in exposures:
            unweighted_RFU.append(temp_HDX_df[temp_HDX_df['log_t']==time]['%d'].to_numpy())
        unweighted_RFU = np.row_stack(unweighted_RFU)

        #weighted_uptake = np.dot(time_weights, unweighted_RFU) / 100
        weighted_uptake = np.mean(unweighted_RFU, axis=0) / 100

        temp_HDX_df = temp_HDX_df.drop_duplicates(subset=['sequence'], keep = 'last')
        start, end = temp_HDX_df['start'].to_numpy(), temp_HDX_df['end'].to_numpy()
        x_label = [f'{chain}_{start[i]}-{end[i]}' for i in range(len(start))]
        return weighted_uptake, x_label
    
    HDX_df = pd.read_excel(HDX_fpath)
    true_apo = {}
    apo_uptake, x_label = [], []

    for protein, state, letter in apo_states:
        uptake, label = get_weighted_uptake(HDX_df, protein, state, letter)
        apo_uptake.append(uptake)
        x_label += label
    apo_uptake = np.concatenate(apo_uptake)

    for i in range(len(x_label)):
        true_apo[x_label[i]] = apo_uptake[i]

    true_complex = {}
    complex_uptake, x_label = [], []

    for protein, state, letter in complex_states:
        uptake, label = get_weighted_uptake(HDX_df, protein, state, letter)
        complex_uptake.append(uptake)
        x_label += label
    complex_uptake = np.concatenate(complex_uptake)

    for i in range(len(x_label)):
        true_complex[x_label[i]] = complex_uptake[i]

    true_diff = {}
    for key in true_apo.keys():
        true_diff[key] = true_complex[key] - true_apo[key]
    return true_diff

def sigmoid(x):
    return 1 / (1 + np.exp(-x))
def binary_cross_entropy(y_true, y_pred):
    epsilon = 1e-12
    y_pred = np.clip(y_pred, epsilon, 1. - epsilon)  # Clip predictions to avoid log(0)
    bce = -np.sum(y_true * np.log(y_pred) + (1-y_true) * np.log(1-y_pred)) / y_true.shape[0]
    return bce
def mean_squared_error(y_true, y_pred, error_limit=0.025):
    return np.mean(((y_true - y_pred)/error_limit) ** 2)
def get_loss(df, batch_name, true_diff, error_limit=0.025):
    batch_df = df[df['Batch'] == batch_name]
    print(f'batch_df: {batch_df.shape}')
    pred_dict = {}
    pred_diff = {}

    for i, row in batch_df.iterrows():
        pred_dict[row['Range']] = row['Y_Pred']

    complex_baseline = df[df['Batch'] == 'finetune'] 
    for i, row in complex_baseline.iterrows(): # FIXME: only for chain A 
        range_label = row['Range']
        if range_label in pred_dict:
            pred_diff[range_label] = pred_dict[range_label] - row['Y_Pred']

    '''
    for chain in batch_df['Chain'].unique(): # only for homomultimer protein, get average for chains
            chain_df = batch_df[batch_df['Chain'] == chain]
            for i, row in chain_df.iterrows():
                range_label = row['Range'].split('_')[1]
                if range_label not in pred_dict.keys():
                    pred_dict[range_label] = []
                pred_dict[range_label].append(float(row['Y_Pred']))
    for key in pred_dict.keys():
        pred_dict[key] = np.mean(pred_dict[key])

    baseline_df = df[df['Batch'] == 'finetune']
    for i, row in baseline_df.iterrows():
        range_label = row['Range'].split('_')[1]
        if range_label in pred_dict:
            pred_diff[f'A_{range_label}'] = pred_dict[range_label] - row['Y_Pred']
    '''

    trim_true_diff = {}
    for key in true_diff.keys():
        if key in pred_diff.keys():
            trim_true_diff[key] = true_diff[key]
    print(f'matched pep for calculating MSE: {len(trim_true_diff.keys())}')
    x_labels = list(trim_true_diff.keys())
    x_index = np.arange(len(x_labels))
    y_true = np.array(list(trim_true_diff.values()))
    y_pred = np.array([pred_diff[key] for key in x_labels])

    loss = mean_squared_error(y_true, y_pred, error_limit)
    return loss

def get_RMSD(pdb_dir, model2, ref_model, ref_chains=['A']):
    parser = PDBParser()
    structure1 = parser.get_structure("Protein1", f"{pdb_dir}/{ref_model}.pdb")
    structure2 = parser.get_structure("Protein2", f"{pdb_dir}/{model2}.pdb")

    # Select the backbone atoms from each structure
    def get_backbone_atoms(structure):
        atoms = []
        for model in structure:
            for chain in model:
                if chain.id in ref_chains:
                    continue
                for residue in chain:
                    try:
                        #atoms.append(residue['N'])
                        atoms.append(residue['CA'].get_coord())
                        #atoms.append(residue['C'])
                    except KeyError:
                        continue
        return atoms

    backbone_atoms1 = get_backbone_atoms(structure1)
    backbone_atoms2 = get_backbone_atoms(structure2)
    # calcualte the distance between backbone_atoms1 and backbone_atoms2
    rmsd = np.sqrt(np.sum((np.array(backbone_atoms1) - np.array(backbone_atoms2))**2)/len(backbone_atoms1))
    # Check if rmsd contains inf or NaN value
    if np.isinf(rmsd) or np.isnan(rmsd):
        print(f"Warning: RMSD for {model2} contains inf or NaN value.")
    return rmsd

def main(config):
    losses = []
    rmsd_list = []
    jaccard_rec= []
    jaccard_lig= []

    # calculate the MSE loss
    batches = config['model_list']
    pdb_dir = config['pdb_dir']
    losses = []
    for batch in batches:
        loss = get_loss(config['prediction'], batch, config['true_diff'], config['error_limit'])
        losses.append(loss)
    #sorted_batch = sorted(zip(batches, losses), key=lambda x: x[1])

    '''# plot the MSE loss
    plt.plot(np.arange(len(batches)), maes)
    plt.xlabel('Model index')
    plt.ylabel('MSE')
    plt.show()
    '''

    '''# get Rosetta score
    Rscore_path = [f'{pdb_dir}/scores/score_{i}.sc' for i in range(1, 101)]
    Rscores = []
    for path in Rscore_path:
        with open(path, 'r') as f:
            lines = f.readlines()
            score = lines[-1].split()[1]
            Rscores.append(float(score))
    
    # get Hdock score
    Hdock_scores = []
    for model in batches:
        pdb_fpath = f'{pdb_dir}/{model}.pdb'
        with open(pdb_fpath, 'r') as f:
            lines = f.readlines()
            score = lines[3].split()[-1]
            Hdock_scores.append(float(score))

    # get structure RMSD and binding site jaccard index
    model_list = config['model_list']
    ref_model = config['ref_model']
    ref_chains = config['ref_chains']

    chains = read_PDB('', f'{pdb_dir}/{ref_model}.pdb')
    _, contact_list1 = bindingsite_extract(chains, chain1='A', chain2='B', dist_cutoff=10.0)

    for model in model_list:
        chains = read_PDB('', f'{pdb_dir}/{model}.pdb')
        _, contact_list2 = bindingsite_extract(chains, chain1='A', chain2='B', dist_cutoff=10.0)
        if contact_list2 is None:
            jaccard_rec.append(0)
            jaccard_lig.append(0)
        else:
            # calculate the jaccard index
            intersection = len(set.intersection(set(contact_list1[:,0]), set(contact_list2[:,0])))
            union = len(set.union(set(contact_list1[:,0]), set(contact_list2[:,0])))
            jaccard = intersection / union
            jaccard_rec.append(jaccard)
            
            intersection = len(set.intersection(set(contact_list1[:,1]), set(contact_list2[:,1])))
            union = len(set.union(set(contact_list1[:,1]), set(contact_list2[:,1])))
            jaccard = intersection / union
            jaccard_lig.append(jaccard)

        rmsd = get_RMSD(pdb_dir, model, ref_model, ref_chains)
        rmsd_list.append(rmsd)
    '''

    '''# plot the RMSD
    plt.figure(figsize = (25, 6))
    plt.plot(np.arange(len(model_list)), rmsd_list)
    plt.xticks(np.arange(len(model_list)), np.arange(len(model_list)))
    plt.grid(True)
    plt.show()
    '''

    '''# calculate the correlation between MSE and RMSD
    jaccard_rec = np.array(jaccard_rec)
    jaccard_lig = np.array(jaccard_lig)
    jaccard_rec = 1 - jaccard_rec
    jaccard_lig = 1 - jaccard_lig

    rescale_RMSD = (rmsd_list - np.min(rmsd_list))/(np.max(rmsd_list) - np.min(rmsd_list))
    rescale_jaccard_rec = (jaccard_rec - np.min(jaccard_rec))/(np.max(jaccard_rec) - np.min(jaccard_rec))
    rescale_jaccard_lig = (jaccard_lig - np.min(jaccard_lig))/(np.max(jaccard_lig) - np.min(jaccard_lig))
    rescale_losses = (losses - np.min(losses))/(np.max(losses) - np.min(losses))
    rescale_Rscores = (Rscores - np.min(Rscores))/(np.max(Rscores) - np.min(Rscores))
    rescale_Hscores = (Hdock_scores - np.min(Hdock_scores))/(np.max(Hdock_scores) - np.min(Hdock_scores))

    combined_struct_diff = rescale_RMSD + rescale_jaccard_rec #+ rescale_jaccard_lig
    #combined_struct_diff = rescale_jaccard_rec + rescale_jaccard_lig
    combined_score = losses
    #combined_struct_diff = np.nan_to_num(combined_struct_diff, nan=2)
    pcc = pearsonr(combined_score, combined_struct_diff)
    scc = spearmanr(combined_score, combined_struct_diff)
    print(pcc)
    print(scc)

    # plot the correlation between MSE and RMSD
    plt.xlabel('MSE between true and pred HDX diff')
    plt.ylabel('RMSD between true and pred structure')
    plt.grid()
    plt.title(f'PCC={pcc[0]:.2f}, SCC={scc[0]:.2f}')
    plt.scatter(losses, rmsd_list)
    plt.show()
    
    return combined_score, combined_struct_diff
    '''

    return losses

if __name__ == '__main__':
    protein_name = 'Hdock_Wuhan+V16ext'

    HDX_fpath = '/home/lwang/models/HDX_LSTM/data/COVID_SPIKE/HDX_files/COVID_SPIKE.xlsx'
    apo_states = [('WUHAN_2nd', 'apo', 'A')] #(protein, state, chain)
    complex_states = [('WUHAN_2nd', 'VH16_VL104', 'A')]
    true_diff = get_true_diff(HDX_fpath, apo_states, complex_states)
    #loss_data = []
    #rmsd_data = []

    #prediction file
    fpath = f'/home/lwang/models/HDX_LSTM/results/240411_BiLSTMGAT_v2_1/HDX_pred_{protein_name}_finetune.csv'
    df = pd.read_csv(fpath)
    config = {
        'prediction': df,
        'true_diff': true_diff,
        'pdb_dir': f'/home/lwang/models/HDX_LSTM/data/COVID_SPIKE/hdock/{protein_name}',
        'model_list': [f'{protein_name}_{i}_revised' for i in range(1, 101)],
        'ref_model': 'wuhan_spikeprotein',
        'ref_chains': ['A'],
        'error_limit': 0.025
    }
    loss_data = main(config)
    print(loss_data)
    model_list = [f'model_{i}' for i in range(1, 101)]
    model_rank = sorted(zip(model_list, loss_data), key=lambda x: x[1])
    for model, loss in model_rank[:10]:
        print(f'{model}: {loss:.2f}')

    ''' # divide the data into 4 groups according to loss value and plot
    loss_bins = [0, 2.2, 2.4, 2.6, 2.8, 4]
    loss_labels = [f'{loss_bins[i]}-{loss_bins[i+1]}' for i in range(len(loss_bins)-1)]
    loss_group = np.digitize(np_loss, bins=loss_bins)

    rmsd_digitize = []
    loss_digitize = []
    for i in range(1, len(loss_bins)):
        rmsd_digitize.append(np_rmsd[loss_group == i])
        loss_digitize.append(np_loss[loss_group == i])
    rmsd_digitize = np.array(rmsd_digitize)
    loss_digitize = np.array(loss_digitize)

    # plot the correlation between MSE and RMSD
    colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'])
    plt.figure()
    for i in range(len(loss_bins)-1):
        #plt.scatter(loss_digitize[i], rmsd_digitize[i], label=f'G{i+1}', color=next(colors))
        plt.scatter(loss_digitize[i], rmsd_digitize[i], color='b')
    plt.xlabel('chi square')
    plt.ylabel('RMSD + Jaccard index')
    plt.title(f'PCC={total_pcc[0]:.2f}, SCC={total_scc[0]:.2f}')
    plt.grid()
    #plt.legend()
    plt.savefig('/home/lwang/models/HDX_LSTM/results/240411_BiLSTMGAT_v2_1/loss_rmsd_correlation_8.png')
    print('scatter Plot saved.')

    # Box plot for loss_data
    plt.figure()
    plt.boxplot(loss_digitize)
    plt.title('Box plot for loss_data')
    plt.xlabel('Index')
    plt.ylabel('RMSD + Jaccard index')
    plt.savefig('/home/lwang/models/HDX_LSTM/results/240411_BiLSTMGAT_v2_1/loss_data_boxplot_8.png')
    print('loss Box plot saved.')

    # Box plot for rmsd_data
    plt.figure()
    plt.boxplot(rmsd_digitize)
    plt.title('Box plot for rmsd_data')
    plt.xlabel('Index')
    plt.ylabel('RMSD + Jaccard index')
    plt.savefig('/home/lwang/models/HDX_LSTM/results/240411_BiLSTMGAT_v2_1/rmsd_data_boxplot_8.png')
    print('rmsd Box plot saved.')
    '''

