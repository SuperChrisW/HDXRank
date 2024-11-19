### structure clustering ###
# 1. clustering based on RMSD

import os
import sys
sys.path.append('/home/lwang/AI-HDX-main/ProteinComplex_HDX_prediction/BiLSTM+GAT_v2/pepGraph_generation')
import numpy as np
import networkx as nx
from pepGraph_dataset import parse_PDB

'''
compute rmsd matrix for HDX score top N decoys, start from no.1 decoys those decoy pairs with rmsd < cutoff will be clustered together 
'''
N=50
rmsd_cutoff = 20.0 # use a coarse cutoff that allows different orientations at a close binding site
atom_type_list = ['CA']
ligand_chain = ['I']

Top_N_dir = '/home/lwang/models/HDX_LSTM/data/hdock/Top_N_hdx'
decoys_coord = []
file_list = []
for file in os.listdir(Top_N_dir):
    if file.startswith(protein_name):
        src_fpath = f'{Top_N_dir}/{file}'
        nodes = parse_PDB(src_fpath, ligand_chain, atom_type_list)
        coord = [res[1]['residue_coord'] for res in nodes]
        coord_array = np.vstack(coord)

        if decoys_coord:
            assert coord_array.shape[0] == decoys_coord[0].shape[0]
        file_list.append(file)
        decoys_coord.append(coord_array)

rmsd_mtx = np.zeros((len(decoys_coord), len(decoys_coord)))
for i in range(len(decoys_coord)):
    for j in range(i, len(decoys_coord)):
        if i == j:
            rmsd_mtx[i, j] = 0.0
        else:
            dist = np.sqrt(np.sum((decoys_coord[i] - decoys_coord[j]) ** 2, axis=1))
            rmsd_mtx[i, j] = rmsd_mtx[j, i] = np.mean(dist)

connect_map = rmsd_mtx < rmsd_cutoff
G = nx.Graph(connect_map)
clusters = list(nx.connected_components(G))
cluster_plot = {}
for cluster_id, cluster in enumerate(clusters):
    print(f'cluster {cluster_id}:')
    cluster_plot[cluster_id] = ([],[])
    mean_HDX_score = 0
    for id in cluster:
        print(file_list[id])
        decoy_id = int(file_list[id][:-4].split('_')[-1])
        mean_HDX_score += HDX_scores[decoy_id-1]
        cluster_plot[cluster_id][0].append(HDX_scores[decoy_id-1])
        cluster_plot[cluster_id][1].append(LRMS[decoy_id-1])
    print('mean HDX score:', mean_HDX_score/len(cluster))


### average translation move and cluster rotations ###

from Bio.PDB import *
from scipy.spatial.transform import Rotation as R 
import numpy as np
import os

def extract_chain_coordinates(pdb_file, chain_id):
    """
    Extracts the coordinates of all atoms in a specified chain from a PDB file.
    
    :param pdb_file: Path to the PDB file
    :param chain_id: ID of the chain to extract
    :return: A NumPy array of shape (n_atoms, 3) containing the atomic coordinates
    """
    structure = parser.get_structure('structure', pdb_file)
    if isinstance(chain_id, str):
        chain_id = [chain_id]

    coordinates = []
    for id in chain_id:
        chain = structure[0][id]
        for atom in chain.get_atoms():
            coordinates.append(atom.coord)
    
    return np.array(coordinates)

def average_coordinates(pdb_files, chain_id, weights=None):
    """
    Averages the coordinates of all atoms in a specified chain across multiple PDB files.
    
    :param pdb_files: List of PDB file paths
    :param chain_id: ID of the chain to average
    :return: A NumPy array of shape (n_atoms, 3) containing the averaged atomic coordinates
    """
    coords_list = []
    
    for pdb_file in pdb_files:
        coords = extract_chain_coordinates(pdb_file, chain_id)
        coords_list.append(coords)
    
    # Stack coordinates into a single array and compute the average
    coords_stack = np.stack(coords_list, axis=0)
    coords_stack = np.mean(coords_stack, axis=1)
    average_coords = np.average(coords_stack, axis=0, weights=weights)
    
    return average_coords
    
# Initialize the parser
parser = PDBParser(QUIET=True)

top_cluster_id = np.argmax([len(cls) for cls in clusters])
cluster_pdb = [file_list[id] for id in clusters[top_cluster_id]]
score = 1-np.array(cluster_plot[top_cluster_id][0])
weights = np.exp(score) / np.sum(np.exp(score)) # softmax, HDX score based weighting / no big difference with uniform weight
print(cluster_pdb)
print(weights)

pdb_fpath = [f'{Top_N_dir}/{file}' for file in cluster_pdb]

lig_chain = ['E']
rec_chain = ['I']
move_factor = 30
rec_center = [average_coordinates([pdb_fpath[0]], chain) for chain in rec_chain]
rec_center = np.mean(rec_center, axis=0)
lig_center = average_coordinates(pdb_fpath, lig_chain, weights=weights)

print(rec_center, lig_center)

moving_lig_center = average_coordinates([pdb_fpath[0]], lig_chain)

sup = Superimposer()
moving_lig_atoms = []
moving_pdb = parser.get_structure('mov', pdb_fpath[0])
for chain in lig_chain:
    moving_lig_atoms.extend(list(moving_pdb[0][chain].get_atoms()))
centered_lig_coords = [atom.get_coord() - moving_lig_center for atom in moving_lig_atoms]

rot_mtx = []
for fpath in pdb_fpath:
    ref_lig_atoms = []
    structure = parser.get_structure('ref', fpath)
    for chain in lig_chain:
        ref_lig_atoms.extend(list(structure[0][chain].get_atoms()))
    sup.set_atoms(ref_lig_atoms, moving_lig_atoms)
    rot_mtx.append(sup.rotran[0])

quaternions = [R.from_matrix(R_mat).as_quat() for R_mat in rot_mtx]
q_avg = np.average(quaternions,axis=0,weights = weights)
q_avg /= np.linalg.norm(q_avg)
R_avg = R.from_quat(q_avg).as_matrix()

superimpose_vec=rec_center-moving_lig_center
apart_direction = rec_center - lig_center
apart_vec = -move_factor * (apart_direction/np.linalg.norm(apart_direction))
trans_mov = superimpose_vec + apart_vec

moving_lig_coords = np.dot(np.array(centered_lig_coords), R_avg) + moving_lig_center + trans_mov

with open(f'{pdb_fpath[0]}', 'r') as f:
    lines = f.readlines()
    with open(f'/home/lwang/models/HDX_LSTM/data/hdock/LocalDock/{protein_name}_localdock.pdb', 'w') as out:
        atom_num = 0
        for line in lines:
            if line.startswith('ATOM'):
                chain_id = line[21]
                if chain_id in lig_chain:
                    line = line[:30] + f'{moving_lig_coords[atom_num][0]:8.3f}{moving_lig_coords[atom_num][1]:8.3f}{moving_lig_coords[atom_num][2]:8.3f}' + line[54:]
                    atom_num += 1
            out.write(line)