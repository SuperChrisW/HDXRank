import os
import copy
import numpy as np
import pandas as pd
import networkx as nx

import torch
import torch.nn.functional as F
from torch.utils.data import Dataset
from torch_geometric.data import Data
from torch_geometric.utils import from_networkx

from Bio.PDB.Polypeptide import protein_letters_3to1
from Bio.PDB import PDBParser, Selection

from pepGraph_utlis import load_embedding, neighbor_search, neighbor_filter


class atoms(object):
    _registry = {}
    def __init__(self, i, name, res_cluster):
        self.__class__._registry[i] = self
        self.i = i
        self.name = name
        self.res_cluster = res_cluster
    
    @classmethod
    def get_atom(cls, atom_index):
        if not isinstance(atom_index, int):
            atom_index = int(atom_index)
        return cls._registry.get(atom_index, None)  # Returns None if atom_index is not found

class res(object):
    _registry = []
    def __init__(self, i, name, chain_id, position):
        self._registry.append(self)
        self.clusters = []
        self.energy = 0
        self.i = i
        self.name = name
        self.chain_id = chain_id
        self.position = position
    
    @classmethod
    def get_res(cls, res_index, chain_id):
        if not isinstance(res_index, int):
            res_index = int(res_index)
        for res in cls._registry:
            if res.i == res_index and res.chain_id == chain_id:
                return res
        return None

class pep(object):
    _registry = []
    #protein_embedding = {}
    def __init__(self, i, pep_sequence, chain, start_pos, end_pos, hdx_value, log_t):
        self._registry.append(self)
        self.i = i
        self.sequence = pep_sequence
        self.chain = chain
        self.start = start_pos
        self.end = end_pos
        self.clusters = []
        self.position = None
        self.node_embedding = None
        self.seq_embedding = None
        self.hdx_value = hdx_value
        self.log_t = log_t

    @classmethod
    def get_pep_adj(cls, linkref, ph_bonds):
        adj = np.zeros([len(cls._registry), len(cls._registry)]) # initialize adjacency matrix
        for so, neighbors in linkref.items():  # linkref record covalent bonds and H bonds
            res1 = atoms.get_atom(so)
            if res1 is None:
                continue
            res1_id = int(res1.res_cluster.split('_')[-1]) # res_cluster => {chain_id}_{res_id}
            res1_chain = res1.res_cluster.split('_')[0]

            pep1 = res.get_res(res1_id, res1_chain)
            pep1_ids = pep1.clusters

            for so_pep in pep1_ids:
                for sf in neighbors:
                    res2 = atoms.get_atom(sf)
                    if res2 is None:
                        continue
                    res2_id = int(res2.res_cluster.split('_')[-1]) # res_cluster => {chain_id}_{res_id}
                    res2_chain = res2.res_cluster.split('_')[0]
                    pep2 = res.get_res(res2_id, res2_chain)
                    pep2_ids = pep2.clusters

                    #if res2_id > len(res._registry):
                    #    continue
                    #pep2_ids = res._registry[res2_id-1].clusters

                    for sf_pep in pep2_ids:
                        if so_pep != sf_pep:
                            adj[so_pep-1, sf_pep-1] = 1
                            adj[sf_pep-1, so_pep-1] = 1
        print('added linkref bonds')

        for phob in ph_bonds: # ph_bonds record hydrophobic interactions
            so = phob[0]
            sf = phob[-1]
            res1 = atoms.get_atom(so)
            if res1 is None:
                continue
            res1_id = int(res1.res_cluster.split('_')[-1]) # res_cluster => {chain_id}_{res_id}
            res1_chain = res1.res_cluster.split('_')[0]

            pep1 = res.get_res(res1_id, res1_chain)
            pep1_ids = pep1.clusters
            for so_pep in pep1_ids:
                res2 = atoms.get_atom(sf)
                if res2 is None:
                    continue
                res2_id = int(res2.res_cluster.split('_')[-1]) # res_cluster => {chain_id}_{res_id}
                res2_chain = res2.res_cluster.split('_')[0]

                pep2 = res.get_res(res2_id, res2_chain)
                pep2_ids = pep2.clusters
                for sf_pep in pep2_ids:
                    if so_pep != sf_pep:
                        adj[so_pep-1, sf_pep-1] = 1
                        adj[sf_pep-1, so_pep-1] = 1
        print('added phob bonds')
        for i in range(len(cls._registry)): # add edge if peptide overlap
            for j in range(i+1, len(cls._registry)):
                if cls._registry[i].end > cls._registry[j].start and cls._registry[i].chain == cls._registry[j].chain:
                    adj[i,j] = 1
                    adj[j,i] = 1
        print('added overlap bonds')
        return torch.Tensor(adj)
    
    @classmethod
    def get_pep_feature(cls, embedding_fpath, max_len = 30): 
        '''
        obtain the feature matrix for peptide
        Note: embedding_fpath is a list of embedding file paths for each chain, {apo_identifier}_{chain}.embedding.txt
        '''
        for fpath in embedding_fpath:
            current_chain = fpath.split('.')[0].split('_')[-1]
            embedding_dict, embedding_matrix, _ = load_embedding(fpath)

            #pep.protein_embedding[current_chain] = embedding_matrix
            print(f'embedding matrix chain {current_chain} shape: {embedding_matrix.shape}')
            for peptide in cls._registry:
                start = peptide.start
                end = peptide.end
                chain = peptide.chain
                if chain == current_chain:
                    mtx = torch.zeros((end-start+1, embedding_matrix.shape[1]), dtype = torch.float32)
                    for index, res_id in enumerate(range(start, end+1)):
                        if res_id in embedding_dict.keys():
                            mtx[index] = embedding_dict[res_id].to(dtype = torch.float32)
                        else:
                            print(f'cannot find embedding for residue {res_id} in chain {chain}')

                    # augment the feat. mtx. with extra features then padding
                    rigidity_std = np.std(mtx[:, 1].view(-1).numpy())
                    seq_length = mtx.shape[0] / max_len
                    extra_vec = (peptide.log_t, rigidity_std, seq_length)

                    extra_col = torch.Tensor(extra_vec).repeat(max_len, 1)
                    padding_needed = max_len - mtx.shape[0]
                    padded_seq = F.pad(mtx, (0, 0, 0, padding_needed), 'constant', 0)
                    extended_padded_seq = torch.cat((padded_seq, extra_col), dim=1)   
                    peptide.seq_embedding = extended_padded_seq

                    mtx_mean = torch.mean(mtx, dim = 0)
                    peptide.node_embedding = torch.cat((mtx_mean, torch.Tensor(extra_vec)), dim = 0)
                    
                    if np.isnan(peptide.node_embedding).any():
                        print(peptide.node_embedding)                        
                        peptide.node_embedding[torch.isnan(peptide.node_embedding)] = 0
                        print(f"Embedding matrix contains NaN values: {fpath}")
                        break

def read_edge(key, proflexdataset_path, H_bond_cutoff = -1.0, central_atom = 'CA'): # read the proflex dataset file and return the bonding information as edge_list
    nhb, nph = -1, -1 # index of hbonds, hydrophobic interactions

    linkref = {} #create linkref list to represent the neighbor atoms/ adj matrix for atoms, including non-covalent bonding
    hbonds = {} # store the informaiton 
    H_energy = {} # store the bond energy
    bond_type = {} # store the non-covalent bond type

    cf_list = [] # central-force bond list
    ph_bonds = [] # hydrophobic interaction , ph_bonds[index] = [[real atom, 3 psedo-atoms index]]
    torsions = [] # record the locked central force pairs
    CA_list = {} # record the CA index
    res._registry = [] # reset the residue registry
    descritpion = {} # store the description of H bond

#    print(path)
    if not os.path.isfile(proflexdataset_path):
        print("cannot find the file", key)
        return None

    with open(proflexdataset_path, 'r') as f:
        data = f.read().strip().split('\n') 
        for line in data:
            #if line[21] is not chain_id: continue # only read the signated chain

            ############## clustering atoms into residue group ###############
            if line[:4] == 'ATOM':
                resn = line[17:20].replace(' ','') # resn is residue name, remove spaces
                natom = int(line[6:11].strip())
                chain_id = line[21]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                position = np.array([x,y,z])
                res_id = line[22:26]

                atoms(natom, resn, f'{chain_id}_{res_id}') #3rd argument is the residue cluster: {chain_id}_{res_id}
                res_id = int(res_id)

                if line[13:15].strip() == central_atom: # should be CA
                    CA_list[res_id] = int(line[6:11])-1 # CA index starts from 0 
                    res(res_id, resn, natom, chain_id, position)


            ############## read covalent bonds except psudoatoms ###############
            elif line[:9] == 'REMARK:CF':
                label, so, sf = line.strip().split()
                so, sf = int(so), int(sf)
                if so <= natom and sf <= natom:
                    cf_list.append([so,sf])

                    if so not in linkref.keys(): # add covalent bonds to the atom edges
                        linkref[so] = []
                    if sf not in linkref.keys():
                        linkref[sf] = []
                    linkref[so].append(sf)
                    linkref[sf].append(so)

            ############## read hydrophobic interaction pairs, ph_bonds[index] = [[real atom, 3 psedo-atoms index]]###############
                elif so <= natom and sf > natom: # atom degree will count after generate ph_bonds list
                    nph += 1
                    ph_bonds.append([so, sf])
                else:
                    ph_bonds[nph].append(sf)

            ############## read locked covalent bonds ###############
            elif line[:9] == 'REMARK:TF': # locked dihedral angles
                if so <= natom and sf <= natom:

                    label, so, sf = line.strip().split()
                    so, sf = int(so), int(sf)
                    torsions.append([sf, so]) # store in a increasing order

            elif line[:9] == 'REMARK:HB':
# label       ID      energy     donor    H    acceptor    description
# REMARK:HB   10     -0.01311    4182    4183    4138    HB Dsp3 Asp2
                nhb += 1
                parts = line.split()
                hb_id, energy, so, s, sf, type = parts[1:7]
                des = " ".join(parts[7:])

                if float(energy) > H_bond_cutoff: continue # only consider the strong H bond

                so, s, sf = int(so), int(s), int(sf)
                if so <= natom and s <= natom and sf <= natom:
                    if nhb not in hbonds.keys():
                        hbonds[nhb] = (so, s, sf)
                        H_energy[nhb] = float(energy)
                        descritpion[nhb] = des
                        if type in ['HB','SB','PH']:
                            bond_type[nhb] = type
                        else:
                            bond_type[nhb] = 'U'

                    # add H bond connection to the atom edges
                    linkref[s].append(sf)
                    linkref[sf].append(s)
                else:
                    ph_bonds[-(nph+1)].append(sf) # [real atom, 3 psedo-atoms, real atom]
                    nph -= 1

            #for bond in ph_bonds:
            #    so = bond[0]
            #    sf = bond[-1]

    return linkref, ph_bonds, torsions, CA_list, natom, (hbonds, H_energy, bond_type), descritpion

def read_HDX_table(HDX_df, proteins, states, chains, correction): # read the HDX table file and return the residue information as node_list
    pep._registry = [] # reset the peptide registry
    npep = 0
    res_chainsplit = {}
    for residue in res._registry:
        if residue.chain_id not in res_chainsplit.keys():
            res_chainsplit[residue.chain_id] = []
        res_chainsplit[residue.chain_id].append(residue)
    print('residue chain split:', res_chainsplit.keys())

    ## seperately process states
    for index, (protein, state, chain) in enumerate(zip(proteins, states, chains)):
        temp_HDX_df = HDX_df[(HDX_df['state']==state) & (HDX_df['protein']==protein)]
        temp_HDX_df = temp_HDX_df.sort_values(by=['start', 'end'], ascending=[True, True])

        print(protein, state, chain, temp_HDX_df.shape)

        for i, row in temp_HDX_df.iterrows():
            sequence = row['sequence'].strip()
            start_pos = int(row['start'])+correction[index]
            end_pos = int(row['end'])+correction[index]
            hdx_value = float(row['%d'])/100
            log_t = row['log_t']

            res_seq = [res for res in res_chainsplit[chain] if start_pos <= res.i <= end_pos]
            pdb_seq = ''.join([protein_letters_3to1[res.name] for res in res_seq])
            if pdb_seq != sequence:
                print("sequence mismatch: chain:", chain, "pdb_Seq:", pdb_seq, 'HDX_seq:', sequence)
                continue
            else:
                peptide = pep(npep, sequence, chain, start_pos, end_pos, hdx_value, log_t) # add chain  
                for residue in res_seq:
                    residue.clusters.append(npep)
                    peptide.clusters.append(residue.i)
                npep += 1
    if npep == 0:
        print("no peptide found")
        return False
    else:
        return True

def distance_ResGraph(pdb_file, distance_cutoff, protein_chain):
    parser = PDBParser()
    structure = parser.get_structure('PDB_structure', pdb_file)
    # Filtering by protein_chain (assuming protein_chain can be a string or a list)
    if isinstance(protein_chain, str):
        chains = [chain for chain in structure.get_chains() if chain.id == protein_chain]
    else:
        chains = [chain for chain in structure.get_chains() if chain.id in protein_chain]
    print('chains:', chains)
    residues = Selection.unfold_entities(chains, 'R')
    G = nx.Graph()
    for i, residue in enumerate(residues):
        residue_id = (residue.get_parent().id, residue.id[1])  # Tuple of (chain ID, residue ID)
        G.add_node(i, residue_id=residue_id)
        res(residue.id[1], residue.get_resname(), residue.get_parent().id, residue['CA'].coord)
    
    # Add edges based on the distance cutoff
    for i, residue1 in enumerate(residues):
        for j, residue2 in enumerate(residues):
            if i >= j:  # Avoid repeating pairs
                continue
            try:
                # Calculate distance between alpha carbons; handle missing 'CA' atoms
                distance = np.linalg.norm(residue1["CA"].coord - residue2["CA"].coord)
                if distance <= distance_cutoff:
                    G.add_edge(i, j, edge_attr=torch.tensor([distance]))
            except KeyError: # Skip residues without a 'CA' atom
                continue
    return G

def rescale(data):
    data_max = np.max(data, axis=0)
    data_min = np.min(data, axis=0)
    
    # Avoid division by zero for constant columns by setting 0/0 to 0 directly
    range = np.where(data_max - data_min == 0, 1, data_max - data_min)
    
    rescaled_data = (data - data_min) / range
    rescaled_data = np.nan_to_num(rescaled_data, nan=0.0)
    
    return rescaled_data

def attach_node_attributes(G, embedding_file):
    embedding_data = torch.load(embedding_file)
    protein_embedding = rescale(embedding_data['embedding'].detach().numpy())

    for i, node in enumerate(G.nodes()):
        embedding = protein_embedding[i]
        G.nodes[node]['x'] = torch.tensor(embedding, dtype = torch.float32)

    return G

class pepGraph(Dataset):
    def __init__(self, keys, root_dir, nfeature, max_len=30, distance_cutoff = 5.0):
        self.keys = keys
        self.root_dir = root_dir
        self.embedding_dir = os.path.join(root_dir, 'embedding_files')
        #self.proflex_dir = os.path.join(root_dir, 'proflex_files')
        self.pdb_dir = os.path.join(root_dir, 'structure')
        self.hdx_dir = os.path.join(root_dir, 'HDX_files')
        self.max_len = max_len
        self.nfeature = nfeature
        self.distance_cutoff = distance_cutoff

    def __len__(self):
        return len(self.keys)

    def __getitem__(self, index):
        '''
        return the graph for the index-th protein complex
        keys format: [database_id, protein, state, match_uni, apo_identifier, chain_identifier, correction, protein_chains]
        '''
        atoms._registry = {}
        res._registry = []
        pep._registry = []

        database_id = self.keys[0][index].strip()
        apo_identifier = self.keys[4][index].strip().split('.')[0]

        protein = self.keys[1][index]
        state = self.keys[2][index]
        match_uni = self.keys[3][index]
        chain_identifier = self.keys[5][index]
        correction = self.keys[6][index]
        protein_chains = self.keys[7][index].split(',')
        print(protein_chains)

        print("processing",database_id, apo_identifier)

        #proflex_fpath = os.path.join(self.proflex_dir, f'{apo_identifier}_clean_Hplus_proflexdataset')
        pdb_fpath = os.path.join(self.pdb_dir, f'{apo_identifier}.pdb')
        target_HDX_fpath = os.path.join(self.hdx_dir, f'{database_id}.xlsx')
        embedding_fpath = os.path.join(self.embedding_dir, f'{apo_identifier}.pt')
        
        #if not os.path.isfile(proflex_fpath):
        #    print("cannot find proflex file: ", proflex_fpath)
        #    return None 
        if not os.path.isfile(target_HDX_fpath):
            print("cannot find HDX table: ", target_HDX_fpath)
            return None 
        if not os.path.isfile(embedding_fpath):
            print("Cannot find embedding file:", embedding_fpath)
            return None
        if not os.path.isfile(pdb_fpath):
            print("Cannot find pdb file:", pdb_fpath)
            return None

        HDX_df = pd.read_excel(target_HDX_fpath, sheet_name='Sheet1')
        ### -------- filter HDX table ------ ####
        HDX_df['state'] = HDX_df['state'].str.replace(' ', '')
        HDX_df['protein'] = HDX_df['protein'].str.replace(' ', '')
        protein = [p.replace(' ', '').strip() for p in protein]
        state = [s.replace(' ', '').strip() for s in state]

        # filter the HDX table
        HDX_df = HDX_df[(HDX_df['state'].isin(state)) & (HDX_df['protein'].isin(protein))]
        HDX_df = HDX_df.drop_duplicates(subset=['sequence', 'state', 'protein'], keep = 'last')
        HDX_df = HDX_df[HDX_df['sequence'].str.len() < self.max_len]
        HDX_df = HDX_df[HDX_df['%d']>0]
        HDX_df = HDX_df.sort_values(by=['start', 'end'], ascending=[True, True]) # mixed states will be separated in read_HDX_table func.
        HDX_df = HDX_df.reset_index(drop=True)
        print('after filtering:', HDX_df.shape)

        # residue graph generation #

        G = distance_ResGraph(pdb_fpath, distance_cutoff=self.distance_cutoff, protein_chain=protein_chains) # generate residue graph
        G = attach_node_attributes(G, embedding_fpath)

        if not read_HDX_table(HDX_df, protein, state, chain_identifier, correction): # read HDX table file
            return None
        print('read in pep numbers:', len(pep._registry))

        def find_neigbhors(nodes, hop=1):
            neigbhor_node = []
            for node in nodes:
                neigbhors = list(G.neighbors(node))
                neigbhor_node.extend(neigbhors)
            if hop > 1:
                neigbhor_node = find_neigbhors(neigbhor_node, hop-1)
            return neigbhor_node
        
        graph_ensemble = []
        i2res_id = {data['residue_id']: node for node, data in G.nodes(data=True)}
        for peptide in pep._registry:
            node_ids = [i2res_id[(peptide.chain, res_id)] for res_id in peptide.clusters if (peptide.chain, res_id) in i2res_id]
            if len(node_ids) == 0:
                continue

            logt_vec = torch.full((len(node_ids), 1), peptide.log_t, dtype=torch.float32)
            seqlen_vec = torch.full((len(node_ids), 1), len(node_ids)/self.max_len, dtype=torch.float32)
            seq_embedding = [G.nodes[node]['x'] for node in node_ids]
            seq_embedding = torch.stack(seq_embedding, dim=0)
            seq_embedding = torch.cat((seq_embedding[:, :self.nfeature], logt_vec, seqlen_vec), dim=1)
            pad_needed = self.max_len - len(node_ids)
            seq_embedding = F.pad(seq_embedding, (0, 0, 0, pad_needed), 'constant', 0)

            neigbhor_node = find_neigbhors(node_ids, hop=1)
            node_ids.extend(neigbhor_node)
            node_ids = set(node_ids)

            subG = G.subgraph(node_ids).copy()
            subG.graph['y'] = peptide.hdx_value
            subG.graph['range'] = (peptide.start, peptide.end)
            subG.graph['chain'] = peptide.chain
            subG.graph['seq_embedding'] = seq_embedding
            data = from_networkx(subG)
            graph_ensemble.append(data)
        label = f'{apo_identifier}'

        return graph_ensemble, label