###########################################################################
##     EXTRACT BINDING SITE FROM PROTEIN-NUCLEIC ACID COMPLEX PDB FILE   ##
##     23/12/15                                                          ##
###########################################################################

from collections import defaultdict
import os
import numpy as np
from scipy.spatial.distance import cdist
from itertools import groupby
import json

class Chain:
    def __init__(self):
        self.atoms = []
        self.sequence_type = None

    def add_atom(self, atom_index, residue_index, residue_type, atom_type, coordinates):
        atom_info = {
            'atom_index': atom_index,
            'residue_index': residue_index,
            'residue_type': residue_type,
            'atom_type': atom_type,
            'coordinates': coordinates
        }
        self.atoms.append(atom_info)
        self._update_sequence_type(residue_type)

    def _update_sequence_type(self, residue_type):
        #amino_acids = {'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 
        #               'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 
        #               'THR', 'VAL', 'TRP', 'TYR'}
        dna_nucleotides = {'DA', 'DT', 'DG', 'DC', 'A', 'T', 'G', 'C'}
        rna_nucleotides = {'DA', 'DU', 'DG', 'DC', 'A', 'U', 'G', 'C'}

        if residue_type in dna_nucleotides:
            self.sequence_type = 'DNA'
        elif residue_type in rna_nucleotides:
            self.sequence_type = 'RNA'
        else:
            self.sequence_type = 'Protein'

    def get_atoms(self):
        return self.atoms

def read_PDB(key, PDB_path):
    if not os.path.isfile(PDB_path):
        print("cannot find the file", key)
        return None
    chains = defaultdict(Chain)

    with open(PDB_path, 'r') as f:
        data = f.read().strip().split('\n') 
        for line in data:
            if line[:4] == 'ATOM':
                n_res = int(line[23:26].strip())
                n_atom = int(line[6:11].strip())
                res_type = line[17:20].strip()
                atom_type = line[12:16].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                chain_id = line[21].strip()
                chains[chain_id].add_atom(n_atom, n_res, res_type, atom_type, [x, y, z])
    return chains

def bindingsite_extract(chains, dist_cutoff=3.65): ###assume only two chains in the PDB file
    ## TODO: given multiple sets of coordinates, calculate the distance matrix
    ChainA_site = []
    ChainB_site = []

    chain_ids = list(chains.keys())

    #print(nuc_chain_id, prot_chain_id)
    coords_a = [atom['coordinates'] for atom in chains[chain_ids[0]].get_atoms()]
    coords_b = [atom['coordinates'] for atom in chains[chain_ids[1]].get_atoms()]

    coords_array_a = np.array(coords_a)
    coords_array_b = np.array(coords_b)

    dist_matrix = cdist(coords_array_a, coords_array_b, 'euclidean')

    pairs = dist_criterion(dist_matrix, cutoff=dist_cutoff) # different criterion can be used here
    for pair in pairs:
        resA_id = pair[0]
        resB_id = pair[1]

        resA = chains[chain_ids[0]].get_atoms()[resA_id]['residue_index']
        resB = chains[chain_ids[1]].get_atoms()[resB_id]['residue_index']
        #nuc_res_type = chains[nuc_chain_id].get_atoms()[nuc_id]['residue_type']
        #prot_res_type = chains[prot_chain_id].get_atoms()[prot_id]['residue_type']

        ChainA_site.append(int(resA))
        ChainB_site.append(int(resB))

    Chain_dict = {}
    Chain_dict[chain_ids[0]] = []
    Chain_dict[chain_ids[0]].extend(sorted(set(ChainA_site)))
    Chain_dict[chain_ids[1]] = []
    Chain_dict[chain_ids[1]].extend(sorted(set(ChainB_site)))

    return Chain_dict

def split_list(a_list, indicator):
    return [list(y) for x, y in groupby(a_list, lambda z: z == indicator) if not x]

def dist_criterion(dist_matrix, cutoff):
    dist_matrix = np.array(dist_matrix)
    mask = dist_matrix < cutoff
    indices = np.where(mask)
    pairs = list(zip(indices[0], indices[1]))
    return pairs

def vectorize(site_list, length):
    ## TODO: generate a [1,0,...] vector of length length, with active site residues as 1 
    vector = np.zeros(length, dtype=int)
    filtered_sites = [site for site in site_list if 0 <= site < length]
    vector[filtered_sites] = 1
    return vector

if __name__ == '__main__':
    PDB_dir = '/Users/liyao/Desktop/Tsuda_Lab/Source_code/Aptamer-binding-prediction/PDB'
    os.chdir(PDB_dir)

    prot_contact_dict = {}
    nuc_contact_dict = {}
    for file in os.listdir(PDB_dir):
        nucleic_dict = {}
        prot_dict = {}
        print(file)
        if file.endswith('.pdb'):
            chains = read_PDB(file, f'{PDB_dir}/{file}')
            prot_site_dict, nuc_site_dict = bindingsite_extract(chains, dist_cutoff=3.65) # 3.65 is the cutoff for protein-nucleic acid contact
            prot_contact_dict[file] = prot_site_dict
            nuc_contact_dict[file] = nuc_site_dict
        #break
    with open('../prot_contact_dict.json', 'w') as f:
        json.dump(prot_contact_dict, f)
    with open('../nuc_contact_dict.json', 'w') as f:
        json.dump(nuc_contact_dict, f)


## OUTPUT FORMAT: PDB_id: {Chain_id: {residue_index, residue_index, ...}, ...}
# prot_contact_dict = {PDB_id: {'A': {1,2,3,4,5,6}, 'B': {1,2,3,4,5,6}}}
# nuc_site_dict = {PDB_id: {'C': {1,2}, 'D': {4,5,6}}}

## NOT FINISHED
# vectorized_dict = { PDB_id: {'A': [1,0,1,0,0,1,0,0,1,0,0,1,0,0,1], 'B': [1,0,1,0,0,1,0,0,1,0,0,1,0,0,1]} }