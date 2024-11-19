import sys
import os
import math
import numpy as np
import pandas as pd

# add amino acid properties matrix to the final vector
## read the matrix table and turn it into a dictionary
### convert sequence to number vector
AA_table = pd.read_csv("/home/lwang/AI-HDX-main/dataset/HDMDvector.csv")
AA_array = AA_table.set_index('ID').T.to_dict('list')
### convert the list in dictionary to array
for key, value in AA_array.items():
    value = list(map(str,value))
    AA_array[key] = value

def MSA_embedding(hhm_file, dssp_file, output_file):
    with open(hhm_file) as f, open(dssp_file) as df, open(output_file ,"w") as out:
        # read the dssp file line by line and add the HMM values in
        sequence_matrix = []
        for line in df:
            hhm_vector = []
            content = line.split()
            hhm_vector = hhm_vector + content[:2] + [content[3]] # extract sequence AA and accessible surface area from dssp.txt
            property_m = AA_array[content[1]] # assign each AA with its own property matrix, adopted from HDMD
            hhm_vector = hhm_vector + property_m # AA property features
            sequence_matrix.append(hhm_vector)
        
        for i in f:
            if i.startswith("#"):
                break
        # start from the lines with HMM values
        for i in range(4):
            f.readline()
        # skip the NULL and HMM headers
        lines = f.read().split("\n")
        # print(len(lines)) ## The file consist of three lines for each AA, first line is the HMM number against each AA,
        ## second line is the 10 conversion values, and the last line is empty. Group the three lines into one AA representative.

        ### skip the - lines i.e. missing residues in fasta files -edited by WANG LIYAO, 2024-09-02
        lines_revised = []
        for idx in range(0,int((len(lines)-2)/3)):
            first_line = lines[idx*3]
            next_line = lines[idx*3+1]
            content1 = first_line.strip().split()
            content2 = next_line.strip().split()
            if content1[0] == '-':
                continue
            else:
                lines_revised.append(first_line)
                lines_revised.append(next_line)
                lines_revised.append('')
        lines_revised.append('')
        lines_revised.append('')
        ### end of skipping - lines

        for idx in range(0,int((len(lines_revised)-2)/3)):
            first_line = lines[idx*3].replace("*","99999") # The * symbol is like NA, so here we assigned a big number to it
            next_line = lines[idx*3+1].replace("*","99999")
            content1 = first_line.strip().split()
            content2 = next_line.strip().split()
            for val1 in content1[2:-1]:
                hhm_val1 = 10/(1 + math.exp(-1 * int(val1)/2000))
                sequence_matrix[idx].append(str(hhm_val1))
            for val2 in content2:
                hhm_val2 = 10/(1 + math.exp(-1 * int(val2)/2000))
                sequence_matrix[idx].append(str(hhm_val2))
        for item in sequence_matrix:
            out.write("\t".join(item)+"\n")

root_dir = '/home/lwang/models/HDX_LSTM/data/hdock'
protein_name = '8F7A'
hhm_dict = {
    'A': 'AF_IMP9_PLDDT70_A',
    'B': 'AF_RAN_PLDDT70_B',
}

hhm_dir = f'{root_dir}/hhm_files'
dssp_dir = f'{root_dir}/dssp_files/{protein_name}'
save_dir = f'{root_dir}/embedding_files/{protein_name}_AIHDX'

if not os.path.exists(save_dir):
    os.makedirs(save_dir)

fail_list = []
for dssp_file in os.listdir(dssp_dir):
    if not dssp_file.endswith('.dssp.txt'):
        continue
    fname = dssp_file.split('.')[0]
    if dssp_file.startswith('MODEL'):
        if fname[-1] in hhm_dict.keys():
            hhm_fpath = os.path.join(hhm_dir, f'{hhm_dict[fname[-1]]}.hhm')
        else:
            fail_list.append(fname)
            continue
    else:
        hhm_fpath = os.path.join(hhm_dir, f'{fname}.hhm')

    dssp_fpath = os.path.join(dssp_dir, dssp_file)
    save_fpath = os.path.join(save_dir, f'{fname}.embedding.txt')
    MSA_embedding(hhm_fpath, dssp_fpath, save_fpath)
print(fail_list)