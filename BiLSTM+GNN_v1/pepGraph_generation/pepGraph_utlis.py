import os
import math
import pandas as pd
import numpy as np
import torch
import torch.nn.functional as F
from pdb2sql import interface
from tool import BioWrappers
from Bio.PDB.Polypeptide import one_to_three as protein_letters_1to3, three_to_one as protein_letters_3to1

import warnings
import pandas as pd
warnings.filterwarnings("ignore")

residue_charge = {'CYS': -0.64, 'HIS': -0.29, 'ASN': -1.22, 'GLN': -1.22, 'SER': -0.80, 'THR': -0.80, 'TYR': -0.80,
                        'TRP': -0.79, 'ALA': -0.37, 'PHE': -0.37, 'GLY': -0.37, 'ILE': -0.37, 'VAL': -0.37, 'MET': -0.37,
                        'PRO': 0.0, 'LEU': -0.37, 'GLU': -1.37, 'ASP': -1.37, 'LYS': -0.36, 'ARG': -1.65}

residue_polarity = {'CYS': 'polar', 'HIS': 'polar', 'ASN': 'polar', 'GLN': 'polar', 'SER': 'polar', 'THR': 'polar', 'TYR': 'polar', 'TRP': 'polar',
                            'ALA': 'apolar', 'PHE': 'apolar', 'GLY': 'apolar', 'ILE': 'apolar', 'VAL': 'apolar', 'MET': 'apolar', 'PRO': 'apolar', 'LEU': 'apolar',
                            'GLU': 'neg_charged', 'ASP': 'neg_charged', 'LYS': 'neg_charged', 'ARG': 'pos_charged'}

polarity_encoding = {'apolar': 0, 'polar': 1, 'neg_charged': 2, 'pos_charged': 3}

ss_list = ['H', 'B', 'E', 'G', 'I', 'T', 'S', 'P', '-']

AA_array = {
    "A": [-0.591, -1.302, -0.733,  1.570, -0.146],
    "C": [-1.343,  0.465, -0.862, -1.020, -0.255],
    "D": [ 1.050,  0.302, -3.656, -0.259, -3.242],
    "E": [ 1.357, -1.453,  1.477,  0.113, -0.837],
    "F": [-1.006, -0.590,  1.891, -0.397,  0.412],
    "G": [-0.384,  1.652,  1.330,  1.045,  2.064],
    "H": [ 0.336, -0.417, -1.673, -1.474, -0.078],
    "I": [-1.239, -0.547,  2.131,  0.393,  0.816],
    "K": [ 1.831, -0.561,  0.533, -0.277,  1.648],
    "L": [-1.019, -0.987, -1.505,  1.266, -0.912],
    "M": [-0.663, -1.524,  2.219, -1.005,  1.212],
    "N": [ 0.945,  0.828,  1.299, -0.169,  0.933],
    "P": [ 0.189,  2.081, -1.628,  0.421, -1.392],
    "Q": [ 0.931, -0.179, -3.005, -0.503, -1.853],
    "R": [ 1.538, -0.055,  1.502,  0.440,  2.897],
    "S": [-0.228,  1.399, -4.760,  0.670, -2.647],
    "T": [-0.032,  0.326,  2.213,  0.908,  1.313],
    "V": [-1.337, -0.279, -0.544,  1.242, -1.262],
    "W": [-0.595,  0.009,  0.672, -2.128, -0.184],
    "Y": [ 0.260,  0.830,  3.097, -0.838,  1.512]
}## adopted from Atchley et al. (2005): https://www.pnas.org/doi/epdf/10.1073/pnas.0408677102

def load_embedding(fpath, feat_rescale = True):
    """
    Load the embedding from a file, return a dictionary of embedding for each residue
    Note: embedding_fpath is a list of embedding file paths for each chain, {apo_identifier}_{chain}.embedding.txt
    embedding file: 0 res_id, 1 res_name, 2-6 HDMD, 7 res_charge, 8-11 res_polarity,
                                    12 SASA, 13 rigidity, 14-16 corrected_hse_mtx, 17- 36 hhm_mtx
    """
    embedding_dict = {}
    embedding = pd.read_table(fpath, header=None, sep=' ')
    embedding_seq = embedding.loc[:,1].to_numpy()
    embedding_index= embedding.loc[:,0].to_numpy()
    res_dict = {id: protein_letters_3to1(res) for id, res in zip(embedding_index, embedding_seq)}

    # sequence-based features
    HDMD = embedding.loc[:, 2:6].to_numpy().reshape(-1, 5)
    Res_charge = embedding.loc[:, 7].to_numpy().reshape(-1, 1)
    Res_polarity = embedding.loc[:, 8:11].to_numpy().reshape(-1, 4)
    ProteinLM = None

    #structure-based features
    SASA = embedding.loc[:, 12].to_numpy().reshape(-1, 1)
    rigidity = embedding.loc[:, 13].to_numpy().reshape(-1, 1)
    HSE = embedding.loc[:, 14:16].to_numpy()

    #MSA-based features
    HHblits = embedding.loc[:, 17:].to_numpy()
    PSSM = None

    #transform
    def rescale(data):
        data_max = np.max(data, axis=0)
        data_min = np.min(data, axis=0)
        rescale_data = (data - data_min[None, :]) / (data_max - data_min)[None, :]
        rescale_data = np.nan_to_num(rescale_data, nan=0.0)
        if feat_rescale:
            return rescale_data
        else:
            return data
    rescale_HSE = rescale(HSE)
    rescale_HHblits = rescale(HHblits)
    rescale_HDMD = rescale(HDMD)
    rescale_Res_charge = rescale(Res_charge)

    # order: SASA, rigidity, HSE, HHblits, HDMD, Res_charge, Res_polarity
    feature_array = np.concatenate((SASA, rigidity, rescale_HSE, rescale_HHblits, 
                                    rescale_HDMD, rescale_Res_charge, Res_polarity), axis=1)

    embedding_dict = {
        embedding.iloc[i, 0]: torch.tensor(feature_array[i, :], dtype = torch.float32) 
        for i in range(embedding.shape[0])}
    embedding_matrix = torch.tensor(feature_array, dtype=torch.float32)

    return embedding_dict, embedding_matrix, res_dict

def seq_embedding(HDX_dataframe, vector_file, pdb_file, protein, state, chain,
                    mismatch_report = True, correction_value = 0, filter = []):
    """
    Generate the input array and truth array for the given HDX data file and embedding file
    filter: list of filters to apply to the data, e.g. contact, unique
    """

    # read HDX data file
    datafile = pd.read_excel(HDX_dataframe, sheet_name='Sheet1')
    datafile.columns = datafile.columns.str.lower()
    datafile = datafile[(datafile['state'].str.strip() == state) & (datafile['protein'].str.strip() == protein)]
    print('after filtering:', datafile.shape)

    # filter peptides > 30 AA
    max_len = 30
    df = datafile[datafile['sequence'].str.len() < max_len] # sequence column
    df = df.sort_values(by=['start', 'end'], ascending=[True, True])

    # keep only the contact residues between chains
    # db -> pdb2sql database
    if 'contact' in filter:
        db = interface(pdb_file)
        cutoff = 8.5
        res_contact_pairs = db.get_contact_residues(
                cutoff=cutoff, return_contact_pairs=True)
        df, keys_to_pop = contact_filter(res_contact_pairs, chain, df)
        print(f"Removed {len(keys_to_pop)} no-contact peptides by distance filter {cutoff}A.")

    if 'unique' in filter:
        df = df.drop_duplicates(subset=['start', 'end'], keep='last')
        print('current shape:', df.shape)

    start_pos = [x + correction_value for x in df['start'].astype(int).tolist()]
    end_pos = [x + correction_value for x in df['end'].astype(int).tolist()]
    seq = df['sequence'].tolist()
    log_t = df['log_t'].to_numpy()

    # read embedding file
    embed_dict, embed_array, res_dict = load_embedding(vector_file, feat_rescale=False)

    row_size = len(start_pos)
    nfactor = embed_array.shape[-1] + 3
    input_array = np.zeros((row_size,max_len,nfactor), dtype = np.float32)

    # filter out the sequences that do not match the embedding file
    pop_idx = []

    for i, (start, end, sequence) in enumerate(zip(start_pos, end_pos, seq)):
        pdb_seq = [res_dict[i] for i in range(start, end+1) if i in res_dict]
        embed_sequence = ''.join(pdb_seq)
        sequence = sequence.replace(' ','')
        if sequence != embed_sequence:
            pop_idx.append(i)
            if mismatch_report:
                print(f"Sequence mismatch at index {i}: {sequence} != {embed_sequence}")
            continue
        else:
            # augment the feat. mtx. with extra features then padding
            seq_array = embed_array[start - 1:end]
            rigidity_std = np.std(seq_array[:, 1].view(-1).numpy())
            seq_length = seq_array.shape[0] / max_len
            extra_vec = (log_t[i], rigidity_std, seq_length)
            extra_column = torch.Tensor(extra_vec).repeat(max_len, 1)

            padding_needed = max_len - seq_array.shape[0]
            padded_seq = F.pad(seq_array, (0, 0, 0, padding_needed), 'constant', 0)
            extended_padded_seq = torch.cat((padded_seq, extra_column), dim=1)            
            input_array[i,:,:] = extended_padded_seq

    start_pos = [x for i, x in enumerate(start_pos) if i not in pop_idx]
    end_pos = [x for i, x in enumerate(end_pos) if i not in pop_idx]
    log_t = [x for i, x in enumerate(log_t) if i not in pop_idx]

    #filter those failed sequences from the input ensembles
    non_empty_mask = ~(input_array == 0).all(axis=(1, 2))
    output_array = input_array[non_empty_mask, :, :] 

    truth = df["%d"].to_numpy()
    truth = truth[non_empty_mask]
    return output_array, truth, start_pos, end_pos, log_t

def contact_filter(res_contact_pairs, chain, df):
    keys_to_pop = []
    for i, row in df.iterrows():
        res_range = range(row['start'], row['end'] + 1)
        if chain == 'A':
            if not any([res[1] in res_range for res in res_contact_pairs.keys()]):
                keys_to_pop.append(i)
        elif chain == 'B':
            contact_list = list(res_contact_pairs.values())
            contact_list = [res for sub_list in contact_list for res in sub_list]
            if not any([res[1] in res_range for res in contact_list]):
                keys_to_pop.append(i)
    df = df.drop(keys_to_pop)
    return df, keys_to_pop

def get_seq_polarity(seq):
    encode_index = [polarity_encoding[residue_polarity[protein_letters_1to3(res)]] for res in seq]
    polarity_mtx = np.zeros((len(seq), 4))
    for i, idx in enumerate(encode_index):
        polarity_mtx[i, idx] = 1
    return polarity_mtx

def one_of_k_encoding(x, allowable_set):
    if x not in allowable_set:
        raise Exception("input {0} not in allowable set{1}:".format(x, allowable_set))
    return [str(int(x == s)) for s in allowable_set]

def load_hhm(hhm_file):
    hhm_mtx = []
    with open(hhm_file) as f:
        for i in f:
            if i.startswith("HMM"):
                break
        # start from the lines with HMM values
        for i in range(2):
            f.readline()
        lines = f.read().split("\n")
        # print(len(lines)) ## The file consist of three lines for each AA, first line is the HMM number against each AA,
        ## second line is the 10 conversion values, and the last line is empty. Group the three lines into one AA representative.
        sequence = ""
        for idx in range(0,int((len(lines)-2)/3)+1):
            first_line = lines[idx*3].replace("*","99999") # The * symbol is like NA, so here we assigned a big number to it
            next_line = lines[idx*3+1].replace("*","99999")
            content1 = first_line.strip().split()
            content2 = next_line.strip().split()
            if content1[0]=='//':
                break
            elif content1[0]=='-':
                continue
            sequence += str(content1[0])
            hhm_val1 = [10/(1 + math.exp(-1 * int(val1)/2000)) for val1 in content1[2:-1]]
            hhm_val2 = [10/(1 + math.exp(-1 * int(val2)/2000)) for val2 in content2]
            hhm_val = hhm_val1 + hhm_val2
            hhm_mtx.append(hhm_val)

    return np.array(hhm_mtx), sequence

def generate_embedding(hhm_file, dssp_file, rigidity_file, pdb_file, output_file, chain_id):
    '''
    Generate embedding file from the pre-computed HMM file, DSSP file, and rigidity file
    '''
    ### convert the list in dictionary to array
    if os.path.isfile(output_file):
        return True
    for key, value in AA_array.items():
        value = list(map(str,value))
        AA_array[key] = value

    rigidity_df = pd.read_table(rigidity_file, sep=' ', header=None)
    rigidity_df.columns = ['residue_index', 'residue_name', 'chain_index', 'rigidity']
    rigidity_df = rigidity_df[rigidity_df['chain_index'] == chain_id]

    dssp_df = pd.read_csv(dssp_file, sep='\t', header=None)
    dssp_df[3] = dssp_df[3].fillna(0.0)
    pdb_seq = dssp_df[1].to_list()
    max_len = dssp_df.shape[0]
    rigidity_df = rigidity_df.iloc[:max_len]
    # make sure the length of rigidity file is the same as dssp file
    # dssp output might lose some residues at the end of the protein
    # dssp file index is not real residue index    
    model = BioWrappers.get_bio_model(pdb_file)  
 
    # sequence label
    res_id = rigidity_df['residue_index'].to_numpy().reshape(-1, 1)
    res_name = rigidity_df['residue_name'].to_numpy().reshape(-1, 1)

    # sequence-based features
    res_charge = np.array([residue_charge[protein_letters_1to3(res)] for res in pdb_seq]).reshape(-1, 1)
    res_polarity = get_seq_polarity(pdb_seq)
    HDMD = np.array([AA_array[res] for res in pdb_seq]).reshape(-1, 5)

    # structure-based features
    hse_dict = BioWrappers.get_hse(model, chain_id)
    SASA = dssp_df[3].to_numpy().reshape(-1, 1)
    rigidity = rigidity_df['rigidity'].to_numpy().reshape(-1, 1)

    # MSA-based features
    hhm_mtx, hhm_seq = load_hhm(hhm_file)

    pdb_seq = ''.join(pdb_seq)
    if hhm_seq == pdb_seq:
        #print('Sequence match between HMM and DSSP')
        pass
    elif hhm_seq[:-1] == pdb_seq:
        hhm_mtx = hhm_mtx[:-1]
    else:
        print(hhm_seq)
        print(pdb_seq)
        raise ValueError('Sequence mismatch between HMM and DSSP')

    corrected_hse_mtx = np.zeros((max_len, 3))
    for i, res_j in enumerate(res_id):
        res_j = str(res_j)
        if (chain_id, res_j) in hse_dict.keys():
            corrected_hse_mtx[i, :] = list(hse_dict[(chain_id, res_j)])
    '''
    print('protein length:', res_id.shape[0])
    print('rigidity length:', rigidity.shape)
    print('dssp length:', SASA.shape)
    print('hse length:', corrected_hse_mtx.shape)
    print('HDMD length:', HDMD.shape)
    print('res_charge length:', res_charge.shape)
    print('res_polarity length:', res_polarity.shape)
    print('hhm length:', hhm_mtx.shape)
    '''
    sequence_mtx = np.concatenate((res_id, res_name, HDMD, res_charge, res_polarity,
                                    SASA, rigidity, corrected_hse_mtx, hhm_mtx), axis=1)
    print('sequence_mtx shape:', sequence_mtx.shape)
    np.savetxt(output_file, sequence_mtx, fmt='%s')

def neighbor_search(node_list, edge_index):
    '''
    Find the neighbors of each node in the graph, return a list of neighbor node indexes
    node_list: list of int
    edge_index: tensor of shape [n, 2]
    '''
    edge_df = pd.DataFrame(edge_index.numpy(), columns=['Node1', 'Node2'])    
    edge_df = edge_df[edge_df['Node1'].isin(node_list)]
    neighbors = edge_df['Node2'].to_list()
    return set(neighbors)

def neighbor_filter(graph, neighbors):
    ''' 
    trim the graph.x, graph.edge_index, graph.edge_attr according to the neighbors list (node indexes)
    '''
    node_to_keep = list(neighbors)
    node_to_keep = sorted(node_to_keep)
    reindex_dict = {node: i for i, node in enumerate(node_to_keep)}

    keep_cols = [graph.x[node].tolist() for node in node_to_keep]
    modified_x = torch.tensor(keep_cols, dtype = torch.float32)

    modified_edge = []
    modified_edge_attr = []
    edge_index = graph.edge_index.t().numpy()
    for index, edge_pair in enumerate(edge_index):
        if edge_pair[0] in node_to_keep and edge_pair[1] in node_to_keep:
            modified_edge.append([reindex_dict[edge_pair[0]], reindex_dict[edge_pair[1]]])
            modified_edge_attr.append(graph.edge_attr[index])
    modified_edge = torch.tensor(modified_edge, dtype=torch.int64)
    modified_edge_attr = torch.tensor(modified_edge_attr, dtype=torch.float32)
    
    graph.x = modified_x
    graph.edge_index = modified_edge.t()
    graph.edge_attr = modified_edge_attr
    return graph