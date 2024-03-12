"""
Created on July 21, 2021
Author: Jiali
Usage: Take the HHBlits output hhm format, and encode the protein sequence
python3 environment
Run:
python MSA_embedding.py <proteinID.hhm> <pdb.dssp.txt> <output txt file for embedding vector>
"""
import os
import math
import numpy as np
import pandas as pd
from Bio import PDB
import warnings
from pepGraph_utlis import generate_embedding

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

ss_list = ['H', 'B', 'E', 'G', 'I', 'T', 'S', 'P', '-']
def one_of_k_encoding(x, allowable_set):
    if x not in allowable_set:
        raise Exception("input {0} not in allowable set{1}:".format(x, allowable_set))
    return [str(int(x == s)) for s in allowable_set]

def generate_embedding(hhm_file, dssp_file, rigidity_file, output_file, chain_id):

    ### convert the list in dictionary to array
    for key, value in AA_array.items():
        value = list(map(str,value))
        AA_array[key] = value

    rigidity_df = pd.read_table(rigidity_file, sep=' ', header=None)
    rigidity_df.columns = ['residue_index', 'residue_name', 'chain_index', 'rigidity']
    rigidity_df = rigidity_df[rigidity_df['chain_index'] == chain_id]
    dssp_df = pd.read_csv(dssp_file, sep='\t', header=None)
    # Replace NAN values in dssp_df column 1 with 0.0
    dssp_df[1] = dssp_df[1].fillna(0.0)    

    max_len = dssp_df.shape[0]
    rigidity_df = rigidity_df.iloc[:max_len]
    correction = max_len - rigidity_df.shape[0]
    # make sure the length of rigidity file is the same as dssp file
    # dssp output might lose some residues at the end of the protein
    # dssp file index is not real residue index

    with open(hhm_file) as f, open(output_file ,"w") as out:
        # read the dssp file line by line and add the HMM values in
        sequence_matrix = []
        for index, row in dssp_df.iterrows():
            hhm_vector = []
            
            dssp_value = round(float(row[3]), 4)
            res_id = rigidity_df.iloc[index]['residue_index'] # use the residue index from rigidity file
            hhm_vector = hhm_vector + [str(res_id), str(row[1]), str(dssp_value)] # extract sequence AA and accessible surface area from dssp.txt
            property_m = AA_array[row[1]] # assign each AA with its own property matrix, adopted from HDMD
            hhm_vector = hhm_vector + property_m # AA property features

            ss_vector = one_of_k_encoding(row[2], ss_list) # secondary structure features
            hhm_vector = hhm_vector + ss_vector

            rigidity = rigidity_df.iloc[index]['rigidity']
            #row = rigidity_df.iloc[index-1]
            rigidity_values = round(float(rigidity),6)
            hhm_vector.append(str(rigidity_values))

            sequence_matrix.append(hhm_vector)

        for i in f:
            if i.startswith("HMM"):
                break
        # start from the lines with HMM values
        for i in range(2):
            f.readline()
        lines = f.read().split("\n")
        # print(len(lines)) ## The file consist of three lines for each AA, first line is the HMM number against each AA,
        ## second line is the 10 conversion values, and the last line is empty. Group the three lines into one AA representative.

        for idx in range(0,int((len(lines)-2)/3)+1):
            first_line = lines[idx*3].replace("*","99999") # The * symbol is like NA, so here we assigned a big number to it
            next_line = lines[idx*3+1].replace("*","99999")
            content1 = first_line.strip().split()
            content2 = next_line.strip().split()

            if idx >= len(sequence_matrix):
                break
            
            for val1 in content1[2:-1]:
                val1 = int(val1)
                hhm_val1 = math.exp(int(val1-5000)/1000)/(1 + math.exp(int(val1-5000)/1000))
                hhm_val1 = round(hhm_val1, 4)
                sequence_matrix[idx].append(str(hhm_val1))
            for val2 in content2:
                val2 = int(val2)
                hhm_val2 = math.exp(int(val2-5000)/1000)/(1 + math.exp(int(val2-5000)/1000))
                hhm_val2 = round(hhm_val2, 4)
                sequence_matrix[idx].append(str(hhm_val2))
            # output the vector for each AA in the protein

        check = np.array(sequence_matrix)
        check = check[:,2:].astype(np.float64)
        print(check.shape)

        if not np.isfinite(check).all():
            print('not finite')

        for item in sequence_matrix:
            out.write("\t".join(item)+"\n")


if __name__ == "__main__":
    ## generate embedding ##
    warnings.filterwarnings("ignore")
    root_dir = '/home/lwang/AI-HDX-main/HDX_MS_dataset/complexpair_dataset'

    dssp_dir = os.path.join(root_dir, 'dssp_files')
    hhm_dir = os.path.join(root_dir, 'hhm_files')
    rigidity_dir = os.path.join(root_dir, 'rigidity_files')
    save_dir = os.path.join(root_dir, 'embedding_files')
    structure_dir = os.path.join(root_dir, 'structrue')
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)

    df = pd.read_excel(f'{root_dir}/merged_complex_pair.xlsx', sheet_name='Sheet1')
    df = df.dropna(subset=['chain_identifier'])
    fail_list = []

    for id, row in df.iterrows():
        apo_pdb = row['apo_identifier'].split('.')[0]
        #if os.path.isfile(f'{save_dir}/{apo_pdb}_{chain_id}.embedding.txt'):
        #    continue
        print('processing:', apo_pdb)

        parser = PDB.PDBParser()
        pdb_file = f'{root_dir}/structure/{apo_pdb}.pdb'
        if not os.path.isfile(pdb_file):
            print('skip')
            continue
        structure = parser.get_structure('pdb', pdb_file)
        chain_ids = list(structure.get_chains())

        for chain_id in chain_ids:
            chain_id = chain_id.get_id()
            dssp_file = os.path.join(dssp_dir, f'{apo_pdb}_{chain_id}.dssp.txt')
            hhm_file = os.path.join(hhm_dir, f'{apo_pdb}_{chain_id}.hhm')
            rigidity_file = os.path.join(rigidity_dir, f'rigidity_{apo_pdb}.txt')

            if os.path.isfile(dssp_file) and os.path.isfile(hhm_file) and os.path.isfile(rigidity_file):
                generate_embedding(hhm_file, dssp_file, rigidity_file, pdb_file,
                                   f'{save_dir}/{apo_pdb}_{chain_id}.embedding.txt', chain_id)
                
            else:
                print(f'file not exist at ', apo_pdb, chain_id)
                fail_list.append((apo_pdb, chain_id))
                continue
