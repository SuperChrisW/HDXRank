import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil

# Comparative modeling by the AutoModel class
from modeller import *              # Load standard Modeller classes
from modeller.automodel.allhmodel import AllHModel
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser
from Bio.Data.PDBData import protein_letters_1to3
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1_extended

import urllib.request
import urllib.error

from Bio import PDB
import csv
import json
import requests
import warnings

def mut_seq(residue_list, Mut_labels):  # ori_seq is a dictonary returned by get_pdb_seq, {residue_idx: residue_name(one letter)}
    # Check if Ori_seq and Mut_labels are provided
    Ori_seq = residue_list.copy()

    if Ori_seq is None or Mut_labels is None:
        print("Error: Ori_seq and Mut_labels are required.")
        return None

    # Check if Mut_labels is a list or array
    if not isinstance(Mut_labels, (list, tuple)):
        print("Error: Mut_labels should be a list or array.")
        return None

    # Iterate over each mutation label
    mutation_made = False  # Flag to track if any mutation has been made

    for Mut_label in Mut_labels:
        # Check if Mut_label is in the correct format
        if len(Mut_label) < 3 or not Mut_label[1:-1].isdigit():
            #print("Error: Invalid Mut_label format.", Mut_label)
            break

        # Extract the mutation position and new amino acid
        position = int(Mut_label[1:-1])
        new_aa = Mut_label[-1]

        # Check if the position is within the sequence length
        if position not in Ori_seq:
            #print(f"Error: Invalid mutation position.", Mut_label)
            break
        else:
            aa = protein_letters_1to3.get(new_aa)
            if aa is not None and is_aa(aa):
                if Ori_seq[position] != Mut_label[0]:
                    #print("Error: the aa mutated cannot match with original sequence: ", Mut_label)
                    break
                else:
                    # Mutate the amino acid at the specified position
                    Ori_seq[position] = new_aa
                    mutation_made = True  # Mutation was made, update the flag
    if not mutation_made:
        return None
    
    # Convert the list back to a string
    mutated_seq = ''
    previous_residue_id = None
    current_residue_id = None
    sorted(Ori_seq.keys())
    for key in Ori_seq.keys():
        aa = protein_letters_1to3.get(Ori_seq[key])
        if aa is not None and is_aa(aa):
            current_residue_id = key
            if previous_residue_id is not None and current_residue_id != previous_residue_id + 1:
                gap_size = current_residue_id - previous_residue_id - 1
                mutated_seq += "-" * gap_size
            previous_residue_id = current_residue_id
            mutated_seq += Ori_seq[key]
#        elif Ori_seq[key] == '.': # '.' represent a non-natural amino acids in HETATOM record
#            current_residue_id = key
#            mutated_seq += Ori_seq[key]
#            previous_residue_id = current_residue_id
        else:
            continue
    return mutated_seq


chain_represent = {'1':'A', '2':'B', '3':'C', '4':'D', '5':'E', '6':'F', '7':'G', '8': 'H', '9':'L', '10':'M', '11': 'N'}

def seqalign_output(seq, chain_num, directory = "/Users/liyao/Desktop/mut_model/seq_alignment/", info = [None], alignment_index = 0):

    #Info = [PDB_wild, MUTATION, Protein_name = None, Source = None, Chain start ID, Chain end ID]
    
    PDB_seq = seq[:chain_num]
    target = seq[-1]
    len_max = len(target)

    PDB_seq_aligned = ""
    target_aligned = ""
    AF_aligned = []

    # Create a PairwiseAligner object
    aligner = PairwiseAligner()

    for sequence in PDB_seq:
    # Perform sequence alignment
        alignment = aligner.align(sequence, target)
        
        # Get the best alignment
        best_alignment = alignment[alignment_index]

        #Splice the seq_record to get the maximum long same sequence, the discontinuous region was covered
        start = int(best_alignment.aligned[1][0][0])
        end = int(best_alignment.aligned[1][-1][-1])
        cut = end - len_max

        AF_aligned.append(best_alignment[1])
        
        target_aligned += "-" * start
        if cut < 0:
            target_aligned += best_alignment[1][start:cut] + "/"
        else:
            target_aligned += best_alignment[1][start:] + "/"

        PDB_seq_aligned += best_alignment[0]+"/"

    PDB_seq_aligned = PDB_seq_aligned[:-1]  #remove the last "/" symbol
    target_aligned = target_aligned[:-1]

    alignment_string = f">P1;"+info[0]+"\n"
    alignment_string += f"structureX:{info[0]}:FIRST:{info[4]}:LAST:{info[5]}::::\n"
    alignment_string += f"{PDB_seq_aligned}*\n\n"
    
    for i in range(chain_num):
        alignment_string += f">P1;"+info[0]+"_AF"+"\n"
        alignment_string += f"structureX:{info[0]}_AF:FIRST:A:LAST:A:{info[3]}:::\n"
        alignment_string += f"{AF_aligned[i]}*\n\n"

    end_chain_letter = chain_represent[str(chain_num)]
    target_string = f">P1;"+info[0]+"_seq"+"\n"
    target_string += f"sequence:{info[0]}_seq:FIRST:A:LAST:{end_chain_letter}:{info[3]}:::\n"
    target_string += f"{target_aligned}*\n"

    alignment_string += target_string

    # Write alignment to a file with .ali suffix
    filename = "alignment_"+info[0]+".ali"
    file_path = os.path.join(directory, filename)

    with open(file_path, "w") as file:
        file.write(alignment_string)
        # print("have created the file")
        return None


chain_represent = {
    '1': 'A', '2': 'B', '3': 'C', '4': 'D', '5': 'E', '6': 'F',
    '7': 'G', '8': 'H', '9': 'I', '10': 'J', '11': 'K', '12': 'L',
    '13': 'M', '14': 'N', '15': 'O', '16': 'P', '17': 'Q', '18': 'R',
    '19': 'S', '20': 'T', '21': 'U', '22': 'V', '23': 'W', '24': 'X',
    '25': 'Y', '26': 'Z'
}

def mutseq_align(PDB_seq, mut_chain_index, directory = "/Users/liyao/Desktop/mut_model/seq_alignment/", info = [None], alignment_index = 0):

    #Info = [PDB_wild, MUTATION, Protein_name = None, Source = None, Chain start ID, Chain end ID]
    chain_num = len(PDB_seq[:-1])
    sequence = PDB_seq[mut_chain_index] #original sequence, ChainA
    target = PDB_seq[-1] # the last one should be mutated seq

    PDB_seq = PDB_seq[:-1]
    len_max = len(target)

    PDB_seq_aligned = ""
    target_aligned = ""

    # Create a PairwiseAligner object
    aligner = PairwiseAligner()

    # Perform sequence alignment
    alignment = aligner.align(sequence, target)
    # Get the best alignment
    best_alignment = alignment[alignment_index]
    #Splice the seq_record to get the maximum long same sequence, the discontinuous region was covered
#    start = int(best_alignment.aligned[1][0][0])
#    end = int(best_alignment.aligned[1][-1][-1])
#    cut = end - len_max
    
    for i in range(chain_num):
        if i == mut_chain_index:
            PDB_seq_aligned += best_alignment[0]+"/"
#                target_aligned += "-" * start
#                if cut < 0:
#                    target_aligned += best_alignment[1][start:cut] + "/"
#                else:
#                    target_aligned += best_alignment[1][start:] + "/"
            target_aligned += best_alignment[1] + "/"
        else:
            target_aligned += PDB_seq[i] + "/"
            PDB_seq_aligned += PDB_seq[i] + "/"
#        print(i)
#        print(PDB_seq_aligned)
    
    PDB_seq_aligned = PDB_seq_aligned[:-1]  #remove the last "/" symbol
    target_aligned = target_aligned[:-1]

    alignment_string = f">P1;{info[0]}\n"
    alignment_string += f"structureX:{info[0]}:FIRST:{info[4]}:LAST:{info[5]}::::\n"
    alignment_string += f"{PDB_seq_aligned}*\n\n"

    end_chain_letter = chain_represent[str(chain_num)]
    alignment_string += f">P1;{info[0]}_{info[1]}\n"
    alignment_string += f"sequence:{info[0]}_{info[1]}:FIRST:A:LAST:{end_chain_letter}::::\n"
    alignment_string += f"{target_aligned}*\n"

    # Write alignment to a file with .ali suffix
    filename = f"alignment_{info[0]}_{info[1]}.ali"
    file_path = os.path.join(directory, filename)

    with open(file_path, "w") as file:
        file.write(alignment_string)
        #print("have created the file")
        return None

'''
protein_letters_1to3 = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe',
    'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu',
    'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
    'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr'
}'''

#protein_letters_3to1["MSE"] = "M"
non_StdAA = ['TRO', 'PYR', 'SNN', 'PTR', 'PCA']
    
def get_pdb_seq(structure, chain_id = 'A'):
    if len(structure) == 0: #if succussfully read the pdb file
        return None, None

    # Get all the chains present in the structure
    chains = list(structure.get_chains())
#    print(f"Number of chains: {len(chains)}")

    for chain in chains:
#    if int(chain_index) < len(chains):
        # Get the specified chain
#        chain_A = chains[chain_index]
#        ID = chain_A.get_id()
#        print(f"Chain ID: {ID}")

        if not chain_id == chain.id:
            continue

        sequence = ""
#        chain_seq = []

        previous_residue_id = None
        residue_indices = {}

        for residue in chain:
            residue_name = residue.get_resname()
            if residue_name == "MSE":
                residue_name = "MET"
#            if residue.id[1] not in residue_indices.keys():
#                residue_indices[residue.id[1]] = []

            if not is_aa(residue_name, standard=False) and residue_name != "MET":
                continue  # Skip the non-amino acid residue
            else:
                if residue_name in non_StdAA:
                    aa = '-'
                else:
                    aa = protein_letters_3to1_extended.get(residue_name)
                idx = int(residue.id[1])
                residue_indices[idx] = protein_letters_3to1_extended.get(residue_name) # this makes sequence have '-' while residue_list transform to normal aa

            # Check for discontinuity
            current_residue_id = residue.get_id()[1]
            if previous_residue_id is not None and current_residue_id != previous_residue_id + 1:
                gap_size = current_residue_id - previous_residue_id - 1
                sequence += "-" * gap_size
            try:
                sequence += aa
            except Exception as e:
                print(residue_name, aa, '//')
                return False, False
            previous_residue_id = current_residue_id
#        chain_seq.append(sequence)

    return sequence, residue_indices

def get_chains(structure):
    if len(structure) != 0: # if successfully read the pdb file
        chains = [chain.id for chain in structure[0].get_chains()]
        return chains
    
def single_model(key, mut_label, dst_dir = '/Users/liyao/Desktop/AF_model', bool_hetatm = False, AF_model = False):

    #################### start single modelling ##############################

    # Change the current working directory
    os.chdir(dst_dir)

    log.verbose()  # Request verbose output
    env = Environ()  # Create a new MODELLER environment to build the model

    # Directory for input atom files
    env.io.atom_files_directory = ['.']

    filename = f"alignment_{key}_{mut_label}"
    env.io.hetatm = bool_hetatm # read in HETATM records from template PDBs

    if AF_model == True:
        knowns = (key, f'{key}_AF')
    else:
        knowns = key
        
    a = AllHModel(env,
            alnfile = filename + '.ali',  # Alignment filename
            knowns = knowns,  # Code of the template
            sequence = f'{key}_{mut_label}')  # Code of the target
    a.starting_model = 1  # Index of the first model
    a.ending_model = 1  # Index of the last model (determines how many models to calculate)

    a.make()  # Perform the actual comparative modeling

def get_AFmodel(uniprot_list, download_dir):
    #download the corresponding model from alphafold database
    os.chdir(download_dir)
    download_list = []
    count = 0
    for key in uniprot_list:
        file = key.split('.')
        key = file[0]
        try:
            url = f"https://alphafold.ebi.ac.uk/files/AF-{key}-F1-model_v4.pdb" ## need uniprot accession
            filename = f"{key}_AF.pdb"
            check_path = os.path.join(download_dir, filename)
            if not os.path.isfile(check_path):
                urllib.request.urlretrieve(url, filename)

            url = f"https://alphafold.ebi.ac.uk/files/AF-{key}-F1-predicted_aligned_error_v4.json" ## error matrix
            filename = f"{key}_AF_error.json"
            check_path = os.path.join(download_dir, filename)
            if not os.path.isfile(check_path):
                urllib.request.urlretrieve(url, filename)
            download_list.append(key)
            count += 1

        except urllib.error.URLError as e:
            print(f"Error occurred while downloading {key} file: {e}")
    print("obtain", count, "PDB files from alphafold database")
    return download_list

def download_pdb_files(pdb_ids, save_dir):
    # URL for downloading PDB files
    url_template = "https://files.rcsb.org/download/{}.pdb"
    
    # Store failed downloads
    failed_downloads = []

    for pdb_id in pdb_ids:
        # Create the URL for this PDB ID
        url = url_template.format(pdb_id)
        print('downloading {}'.format(pdb_id))
        # Send a GET request to the URL
        response = requests.get(url)

        # Check if the request was successful
        if response.status_code == 200:
            # Write the response content to a file
            fpath = os.path.join(save_dir, f"{pdb_id}.pdb")
            with open(fpath, 'wb') as f:
                f.write(response.content)
        else:
            print(f"Failed to download {pdb_id}. Status code: {response.status_code}")
            failed_downloads.append(pdb_id)

    # Write failed downloads to a file
    fpath = os.path.join(save_dir, "fail_download.txt")
    with open(fpath, 'w') as f:
        for pdb_id in failed_downloads:
            f.write(f"{pdb_id}\n")

def correct_thermomut_dataset(AF_dict, wild_dict, update_index_list, dataset, save_fpath):
    warnings.filterwarnings("ignore")
    df = dataset
    rest_list = []

    for idx in update_index_list:
        row_update = df.loc[idx]
        #print(row_update)
    
        PDB = row_update['PDB Id'].upper()
        #print('row index:', idx)
        UNI_id = row_update['UniProt']
        mut_chain = row_update['Mutated Chain']
        UNI_mut_label = row_update['Mutation_UNP']
        if pd.isna(row_update['Mutation_PDB']):
            mut_label = row_update['Mutation_UNP']
        else:
            mut_label = row_update['Mutation_PDB']

        if not pd.isna(UNI_mut_label) and not pd.isna(mut_label): 
            UNI_label_list = [label.strip() for label in UNI_mut_label.split(',')]
            label_list = [label.strip() for label in mut_label.split(',')]
        else: 
            #print(1)
            rest_list.append(idx)
            #break
            continue
        if UNI_id in AF_dict.keys(): 
            AF_seq, AF_seq_list = AF_dict[UNI_id]['A']
            try:
                WT_seq, WT_seq_list = wild_dict[PDB][mut_chain]
            except Exception as e:
                print(PDB, idx, mut_chain)
                continue
            #print(WT_seq_list)
        else:
            #print(2)
            rest_list.append(idx)
            #break
            continue

        label = UNI_label_list[0] #only check the first mut_label
        position = int(label[1:-1])
        ori_aa = label[0]
        new_aa = label[-1]
        if position <= len(AF_seq):
            if not AF_seq_list[position] == ori_aa:
                rest_list.append(idx)
                #print(3)
                #break
                continue
            else:
                aligner = PairwiseAligner()
                alignment = aligner.align(WT_seq, AF_seq)
                try:
                    aligned_wt = alignment[0][0]
                    aligned_af = alignment[0][1]
                    #print(alignment[0])
                    index = get_wt_index_from_af_position(aligned_wt, aligned_af, position)
                    update = int(label_list[0][1:-1])
                    base_index = int(list(WT_seq_list.keys())[0])
                    for i, item in enumerate(label_list):
                        update = int(item[1:-1]) - update
                        updated_item = item[0] + str(index + base_index + update - 1) + item[-1]
                        label_list[i] = updated_item
                    string = ','.join(label_list)
                    
                    row_update['Mutation_PDB'] = string
                    df.loc[idx] = row_update
                    #print(4)
                except Exception as e:
                    print('position:', position)
                    print(UNI_label_list)
                    print(label_list)
                    print(AF_seq_list)
                    print(e)
            #break
        else:
            rest_list.append(idx)
    
    df.to_csv(save_fpath, index = False)
    print('rest: ', rest_list)
    return rest_list


def get_wt_index_from_af_position(aligned_wt, aligned_af, af_position):
    wt_index = 0
    af_index = 0
    
    for i in range(len(aligned_wt)):
        if aligned_wt[i] != '-':
            wt_index += 1
        if aligned_af[i] != '-':
            af_index += 1
            
        # Check if we've reached the desired position in the AF sequence
        if af_index == af_position:  # +1 because indices are 0-based but positions are 1-based
            return wt_index  #1-based index 
    
    return None