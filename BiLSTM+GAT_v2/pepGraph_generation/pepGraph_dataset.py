import os
import numpy as np
import pandas as pd
import networkx as nx

import torch
import torch.nn.functional as F
from torch.utils.data import Dataset
import torch_geometric
from torch_geometric.data import HeteroData
import torch_cluster
from torch_cluster import knn_graph, radius_graph

from torchdrug import data
from torchdrug.layers import geometry

from Bio.PDB.Polypeptide import protein_letters_3to1
from Bio.PDB import PDBParser, Selection
from Bio import BiopythonWarning
import warnings
warnings.filterwarnings("ignore", category=BiopythonWarning)

from GVP_graphdataset import ProteinGraphDataset, _rbf, _normalize

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
    def __init__(self, i, pep_sequence, chain, start_pos, end_pos, hdx_value):
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

def read_HDX_table(HDX_df, proteins, states, chains, correction, protein_chains, cluster_index): # read the HDX table file and return the residue information as node_list
    pep._registry = [] # reset the peptide registry
    npep = 0
    res_chainsplit = {}
    timepoints = [1.35, 2.85] #predefined log_t boundaries according to k-mean clustering results

    for residue in res._registry:
        if residue.chain_id not in res_chainsplit.keys():
            res_chainsplit[residue.chain_id] = []
        res_chainsplit[residue.chain_id].append(residue)
    print('residue chain split:', res_chainsplit.keys())

    ## seperately process states
    for index, (protein, state, chain) in enumerate(zip(proteins, states, chains)):
        temp_HDX_df = HDX_df[(HDX_df['state']==state) & (HDX_df['protein']==protein)]
        temp_HDX_df = temp_HDX_df.sort_values(by=['start', 'end'], ascending=[True, True])
        chain_index = protein_chains.index(chain)
        print(protein, state, chain_index, temp_HDX_df.shape)

        clusters = [temp_HDX_df[temp_HDX_df['log_t'] < timepoints[0]],
                    temp_HDX_df[(temp_HDX_df['log_t'] >= timepoints[0]) & (temp_HDX_df['log_t'] < timepoints[1])],
                    temp_HDX_df[temp_HDX_df['log_t'] >= timepoints[1]]]

        HDX_grouped = temp_HDX_df.groupby(['start', 'end']) 
        for name, group in HDX_grouped:
            sequence = group['sequence'].iloc[0].strip()
            start_pos = int(name[0])+correction[index]
            end_pos = int(name[1])+correction[index]

            cluster_rfu_means = []
            cluster = clusters[cluster_index]
            cluster_group = cluster[(cluster['start']==name[0]) & (cluster['end']==name[1])]
            if not cluster_group.empty:
                mean_rfu = cluster_group['RFU'].mean()/100
            else:
                continue
            cluster_rfu_means.append(mean_rfu)

            res_seq = [res for res in res_chainsplit[chain_index] if start_pos <= res.i <= end_pos]
            pdb_seq = ''.join([protein_letters_3to1[res.name] for res in res_seq])
            if pdb_seq != sequence: # residue lack/mutation/mismatch in pdb will be filtered out
                print("sequence mismatch: chain:", chain, "pdb_Seq:", pdb_seq, 'HDX_seq:', sequence)
                continue
            else:
                peptide = pep(npep, sequence, chain_index, start_pos, end_pos, cluster_rfu_means[0])
                for residue in res_seq:
                    residue.clusters.append(npep)
                    peptide.clusters.append(residue.i)
                npep += 1
    if npep == 0:
        print("no peptide found")
        return False
    else:
        return True

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
    for node in G.nodes():
        G.nodes[node]['x'] = torch.tensor(protein_embedding[node, :], dtype=torch.float32) # revise the feat attached here
    return G

def ResGraph(pdb_file, embedding_file, protein_chain = ['A']): # create networkx and torchdrug graphs
    parser = PDBParser()
    structure = parser.get_structure('PDB_structure', pdb_file)
    # Filtering by protein_chain (assuming protein_chain can be a string or a list)
    if isinstance(protein_chain, str):
        chains = [chain for chain in structure.get_chains() if chain.id == protein_chain]
    else:
        chains = [chain for chain in structure.get_chains() if chain.id in protein_chain]
    residues = Selection.unfold_entities(chains, 'R')

    G = nx.MultiGraph()
    for i, residue in enumerate(residues):
        nodes_with_attributes = [
            (i, {"res_name": residue.get_resname(), "residue_coord":residue['CA'].coord, 
                 "residue_id": residue.id[1], 'chain_id': protein_chain.index(residue.get_parent().id)})
        ]
        G.add_nodes_from(nodes_with_attributes)
        res(residue.id[1], residue.get_resname(), protein_chain.index(residue.get_parent().id), residue['CA'].coord)

    G = attach_node_attributes(G, embedding_file)
    return G

def GVP_ResGraph(pdb_file, embedding_file, protein_chain = ['A']): # create networkx and torchdrug graphs
    parser = PDBParser()
    structure = parser.get_structure('PDB_structure', pdb_file)
    # Filtering by protein_chain (assuming protein_chain can be a string or a list)
    if isinstance(protein_chain, str):
        chains = [chain for chain in structure.get_chains() if chain.id == protein_chain]
    else:
        chains = [chain for chain in structure.get_chains() if chain.id in protein_chain]
    residues = Selection.unfold_entities(chains, 'R')

    G = nx.MultiGraph()
    for i, residue in enumerate(residues):
        Ca_coord = residue['CA'].coord
        C_coord = residue['C'].coord
        N_coord = residue['N'].coord
        O_coord = residue['O'].coord if 'O' in residue else [-1, -1, -1]
        backbone_coord = np.array([N_coord, Ca_coord, C_coord, O_coord])

        nodes_with_attributes = [
            (i, {"res_name": residue.get_resname(), "residue_coord":backbone_coord, 
                 "residue_id": residue.id[1], 'chain_id': protein_chain.index(residue.get_parent().id)})
        ]
        G.add_nodes_from(nodes_with_attributes)

        res(residue.id[1], residue.get_resname(), protein_chain.index(residue.get_parent().id), backbone_coord)
    G = attach_node_attributes(G, embedding_file)
    G.graph['name'] = pdb_file

    return G

def add_edges(G, coord = None, max_distance=10, min_distance=5):
    batch = torch.tensor([0] * len(G.nodes))

    if coord is not None:
        node_coord = torch.as_tensor(coord, dtype=torch.float32)
    else:
        node_coord = np.array([G.nodes[i]['residue_coord'] for i in G.nodes])
        node_coord = torch.tensor(node_coord, dtype = torch.float32)

    # add radius edges
    edge_index = radius_graph(node_coord, r=max_distance, batch=batch, loop=False)
    edge_list = edge_index.t().tolist()
    for u, v in edge_list:
        if u < v:
            G.add_edge(u, v, edge_type='radius_edge')
    print('add radius edges:', G.number_of_edges())

    # add knn edges
    edge_index = knn_graph(node_coord, k=10, batch=batch, loop=False)
    edge_list = edge_index.t().tolist()
    for u, v in edge_list:
        if u < v:
            G.add_edge(u, v, edge_type='knn_edge')
    print('add knn edges:', G.number_of_edges())

    # remove edges with distance smaller than min_distance
    edges_to_remove = []
    for u, v in G.edges():
        u_coord = node_coord[u]
        v_coord = node_coord[v]
        distance = torch.norm(u_coord - v_coord, dim=0)
        if distance < min_distance:
            if (u,v) not in edges_to_remove:
                edges_to_remove.append((u, v))
    G.remove_edges_from(edges_to_remove)
    print('after removing edges:', G.number_of_edges())

    # add sequential edges
    i2res_id = {(data['chain_id'], data['residue_id']): node for node, data in G.nodes(data=True)}
    for node in G.nodes:
        res_id = G.nodes[node]['residue_id']
        chain = G.nodes[node]['chain_id']
        if (chain,res_id+1) in i2res_id.keys():
            G.add_edge(node, node+1, edge_type='forward_1_edge')
        if (chain,res_id+2) in i2res_id.keys():
            G.add_edge(node, node+2, edge_type='forward_2_edge')
        if (chain,res_id-1) in i2res_id.keys():
            G.add_edge(node, node-1, edge_type='backward_1_edge')
        if (chain,res_id-2) in i2res_id.keys():
            G.add_edge(node, node-2, edge_type='backward_2_edge')  
        G.add_edge(node, node, edge_type='self_edge')
    print('add sequential edges:', G.number_of_edges())

    return G

def networkx_to_HeteroG(subG): # convert to pytorch geometric HeteroData
    data = HeteroData()

    id_map = {}
    for i, (node, node_attr) in enumerate(subG.nodes(data=True)):
        # Assuming node feature vector 'x' is already a tensor or can be converted as such
        #data['residue'].x = torch.cat([data['residue'].x, node_attr['x'].unsqueeze(0)], dim=0) if 'x' in data['residue'] else node_attr['x'].unsqueeze(0)
        if node not in id_map.keys():
            id_map[node] = i

    edge_index_dict = {}
    for u, v, edge_attr in subG.edges(data=True):
        edge_type = edge_attr['edge_type']
        edge_label = ('residue', edge_type, 'residue')
        edge_index = [id_map[u], id_map[v]]
        if edge_label not in edge_index_dict.keys():
            edge_index_dict[edge_label] = edge_index
        else:
            edge_index_dict[edge_label].extend(edge_index)
    for edge_label, edge_index in edge_index_dict.items():
        edge_index = torch.as_tensor(edge_index, dtype=torch.long).t().contiguous()
        data[edge_label].edge_index = edge_index

    if 'y' in subG.graph:
        data['residue'].y = torch.as_tensor([subG.graph['y']], dtype=torch.float32)
    if 'range' in subG.graph:
        data['residue'].range = torch.as_tensor([subG.graph['range']], dtype=torch.float32)
    if 'chain' in subG.graph:
        data['residue'].chain = torch.as_tensor([subG.graph['chain']], dtype=torch.int64)
    if 'seq_embedding' in subG.graph:
        data['residue'].seq_embedding = torch.as_tensor(subG.graph['seq_embedding'], dtype=torch.float32)
    if 'is_complex' in subG.graph:
        data['residue'].is_complex = torch.as_tensor([subG.graph['is_complex']], dtype=torch.int64)

    return data

def networkx_to_tgG(G): # convert to torchdrug protein graph
    node_position = torch.as_tensor(np.array([G.nodes[node]['residue_coord'] for node in G.nodes()]), dtype=torch.float32)
    num_atom = G.number_of_nodes()
    atom_type = ['CA'] * num_atom # CA - 1
    atom_type = torch.as_tensor([data.Protein.atom_name2id.get(atom, -1) for atom in atom_type])

    residue_type = []
    residue_feature = []
    residue_id = []
    id_map = {}
    for i, (node, attrs) in enumerate(G.nodes(data=True)):
        if node not in id_map.keys():
            id_map[node] = i
        residue_type.append(data.Protein.residue2id.get(attrs['res_name'], 0))
        residue_feature.append(attrs['x'])
        residue_id.append(attrs['residue_id'])
    residue_type = torch.tensor(residue_type, dtype=torch.long)
    residue_feature = torch.stack(residue_feature)
    atom2residue = torch.as_tensor([id_map[node] for node in G.nodes()], dtype=torch.long)

    edge_list = []
    bond_type = []
    edge_type_list = ['knn_edge', 'radius_edge', 'self_edge', 'forward_1_edge', 'forward_2_edge', 'backward_1_edge', 'backward_2_edge'] #radius_edge
    #'self_edge', 'forward_1_edge', 'forward_2_edge', 'backward_1_edge', 'backward_2_edge'

    for u, v, attrs in G.edges(data=True):
        edge_type = attrs['edge_type']
        u = id_map[u]
        v = id_map[v]
        if edge_type in edge_type_list:
            edge_list.append([u, v, edge_type_list.index(edge_type)])
            bond_type.append(edge_type_list.index(edge_type))

    edge_list = torch.tensor(edge_list, dtype=torch.long)
    bond_type = torch.tensor(bond_type, dtype=torch.long).unsqueeze(-1)


    protein = data.Protein(edge_list, atom_type, bond_type, view='residue', residue_number=residue_id,
                           node_position=node_position, atom2residue=atom2residue,residue_feature=residue_feature, 
                           residue_type=residue_type, num_relation = len(edge_type_list))

    Constructor = geometry.GraphConstruction(edge_feature = 'gearnet') 
    edge_feature = Constructor.edge_gearnet(protein, protein.edge_list, len(edge_type_list))

    with protein.graph():
        protein.y = torch.as_tensor(G.graph['y'], dtype=torch.float32)
        protein.range = torch.as_tensor(G.graph['range'], dtype=torch.float32)
        protein.chain = torch.as_tensor(G.graph['chain'], dtype=torch.int64)
        protein.is_complex = torch.as_tensor(G.graph['is_complex'], dtype=torch.int64)
        protein.seq_embedding = torch.as_tensor(G.graph['seq_embedding'], dtype=torch.float32)
        #protein.pdb = G.graph['pdb']
        #protein.sequence = G.graph['sequence']
    
    with protein.edge():
        protein.edge_feature = edge_feature

    return protein

class modified_GVPdataset(ProteinGraphDataset):
    def _featurize_as_graph(self, protein):
        # from GVP 
        # modify the code to accept networkx graph as input

        name = protein.graph['name']
        with torch.no_grad():
            residue_coord = np.array([protein.nodes[node]['residue_coord'] for node in protein.nodes()])
            coords = torch.as_tensor(residue_coord, dtype=torch.float32)
            x = torch.stack([protein.nodes[node]['x'] for node in protein.nodes()], dim=0)
            mask = torch.isfinite(coords.sum(dim=(1,2)))
            coords[~mask] = np.inf
            
            X_ca = coords[:, 1]
            edges = list(protein.edges())
            edge_index = torch.tensor(edges, dtype=torch.long, device=self.device)
            edge_index = edge_index.t().contiguous()

            pos_embeddings = self._positional_embeddings(edge_index)
            E_vectors = X_ca[edge_index[0]] - X_ca[edge_index[1]]
            rbf = _rbf(E_vectors.norm(dim=-1), D_count=self.num_rbf, device=self.device)
            
            #dihedrals = self._dihedrals(coords)                     
            orientations = self._orientations(X_ca)
            sidechains = self._sidechains(coords)
            
            node_s = x
            node_v = torch.cat([orientations, sidechains.unsqueeze(-2)], dim=-2)
            edge_s = torch.cat([rbf, pos_embeddings], dim=-1)
            edge_v = _normalize(E_vectors).unsqueeze(-2)
            
            node_s, node_v, edge_s, edge_v = map(torch.nan_to_num,
                    (node_s, node_v, edge_s, edge_v))

            for idx, node in enumerate(protein.nodes()):
                protein.nodes[node]['node_s'] = node_s[idx]
                protein.nodes[node]['node_v'] = node_v[idx]

            for idx, (u, v, k) in enumerate(protein.edges(keys=True)):
                protein.edges[(u,v,k)]['edge_s'] = edge_s[idx]
                protein.edges[(u,v,k)]['edge_v'] = edge_v[idx]

        return protein

class pepGraph(Dataset):
    def __init__(self, keys, root_dir, cluster_index, nfeature, max_len=30, min_distance = 5.0 , max_distance = 10.0,
                 protein_name = None, truncation_window_size = None, graph_type = 'GearNet'):
        self.keys = keys
        self.root_dir = root_dir
        self.cluster_index = cluster_index

        hdock_protein = protein_name
        self.graph_type = graph_type

        self.embedding_dir = os.path.join(root_dir, 'embedding_files')
        self.pdb_dir = os.path.join(root_dir, 'structure')
        self.hdx_dir = os.path.join(root_dir, 'HDX_files')
        self.save_dir = os.path.join(root_dir, 'graph_ensemble_GearNetEdge', f'cluster{self.cluster_index}')

        self.max_len = max_len
        self.nfeature = nfeature
        self.min_distance = min_distance
        self.max_distance = max_distance
        self.truncation_window_size = truncation_window_size

        self.complex_state = {'single':0, 'protein complex':1, 'ligand complex':2, 'dna complex':3, 'rna complex': 4}

    def __len__(self):
        return len(self.keys)

    def __getitem__(self, index):
        if self.graph_type == 'GearNet':
            return self._getitem_GearNet(index)
        elif self.graph_type == 'GVP':
            return self._getitem_GVP(index)
        else:
            raise ValueError(f'Unsupported graph type: {self.graph_type}')

    def _getitem_GearNet(self, index):
        '''
        return the graph for the index-th protein complex
        keys format: [database_id, protein, state, match_uni, apo_identifier, chain_identifier, correction, protein_chains]
        '''
        atoms._registry = {}
        res._registry = []
        pep._registry = []

        database_id = self.keys[0][index].strip()
        apo_identifier = self.keys[4][index].strip().split('.')[0].upper()

        protein = self.keys[1][index]
        state = self.keys[2][index]
        match_uni = self.keys[3][index]
        chain_identifier = self.keys[5][index]
        correction = self.keys[6][index]
        protein_chains = self.keys[7][index].split(',')
        complex_state = self.keys[8][index]

        print("processing",database_id, apo_identifier)
        #if os.path.isfile(os.path.join(self.save_dir, f'{apo_identifier}.pt')):
        #    print('skip')
        #    return None

        pdb_fpath = os.path.join(self.pdb_dir, f'{apo_identifier}.pdb')
        target_HDX_fpath = os.path.join(self.hdx_dir, f'{database_id}_revised.xlsx')
        embedding_fpath = os.path.join(self.embedding_dir, f'{apo_identifier}.pt')
        
        if not os.path.isfile(target_HDX_fpath) and self.truncation_window_size is None:
            print("cannot find HDX table: ", target_HDX_fpath)
            return None 
        if not os.path.isfile(embedding_fpath):
            print("Cannot find embedding file:", embedding_fpath)
            return None
        if not os.path.isfile(pdb_fpath):
            print("Cannot find pdb file:", pdb_fpath)
            return None

        # residue graph generation #
        G = ResGraph(pdb_fpath, embedding_fpath, protein_chain=protein_chains) # generate residue graph
        G = add_edges(G, max_distance=self.max_distance, min_distance=self.min_distance)

        if self.truncation_window_size is None:
            HDX_df = pd.read_excel(target_HDX_fpath, sheet_name='Sheet1')
            ### -------- filter HDX table ------ ####
            HDX_df['state'] = HDX_df['state'].str.replace(' ', '')
            HDX_df['protein'] = HDX_df['protein'].str.replace(' ', '')
            protein = [p.replace(' ', '').strip() for p in protein]
            state = [s.replace(' ', '').strip() for s in state]

            # filter the HDX table
            HDX_df = HDX_df[(HDX_df['state'].isin(state)) & (HDX_df['protein'].isin(protein))]
            HDX_df = HDX_df[HDX_df['sequence'].str.len() < self.max_len]
            HDX_df = HDX_df.sort_values(by=['start', 'end'], ascending=[True, True]) # mixed states will be separated in read_HDX_table func.
            HDX_df = HDX_df.reset_index(drop=True)
            print('after filtering:', HDX_df.shape)
            
            if not read_HDX_table(HDX_df, protein, state, chain_identifier, correction, protein_chains, self.cluster_index): # read HDX table file
                return None
            print('read in pep numbers:', len(pep._registry))
        else: # instead of reading from HDX source file, move a truncation window along sequence
            npep = 0
            stride = 1
            chain_batch = []
            res_name = []
            res_id = []
            for residue in res._registry:
                res_id.append(residue.i)
                res_name.append(residue.name)
                chain_batch.append(residue.chain_id)
            res_id, res_name, chain_batch = np.array(res_id), np.array(res_name), np.array(chain_batch)

            for chain_label in set(chain_batch):
                mask = (chain_batch == chain_label)
                chain_res_id = res_id[mask]
                chain_res_name = res_name[mask]
                for i in range(0, len(chain_res_id) - self.truncation_window_size, stride):
                    seq = ''.join([protein_letters_3to1[res] for res in chain_res_name[i:i+self.truncation_window_size]])
                    start_pos, end_pos = chain_res_id[i], chain_res_id[i+self.truncation_window_size-1]
                    hdx_value = -1
                    peptide = pep(npep, seq, chain_label, start_pos, end_pos, hdx_value)
                    for j in range(len(seq)):
                        peptide.clusters.append(chain_res_id[i+j])
                    npep += 1
            print('read in pep numbers:', len(pep._registry))

        def find_neigbhors(G, nodes, hop=1):
            neigbhor_node = []
            for node in nodes:
                for key, neighbor, edge_data in G.edges(node, data=True):
                    if edge_data['edge_type'] == 'radius_edge':
                        neigbhor_node.append(neighbor)
            if hop > 1:
                neigbhor_node = find_neigbhors(neigbhor_node, hop-1)
            return set(neigbhor_node)
        
        i2res_id = {(data['chain_id'], data['residue_id']): node for node, data in G.nodes(data=True)}
        graph_ensemble = []

        for peptide in pep._registry:
            chain = peptide.chain
            node_ids = [i2res_id[(chain, res_id)] for res_id in peptide.clusters if (chain, res_id) in i2res_id]
            if len(node_ids) == 0:
                continue
            
            seq_embedding = [G.nodes[node]['x'] for node in node_ids]
            seq_embedding = torch.stack(seq_embedding, dim=0)
            pad_needed = self.max_len - len(node_ids)
            seq_embedding = F.pad(seq_embedding, (0, 0, 0, pad_needed), 'constant', 0)
            
            neigbhor_node = find_neigbhors(G, node_ids, hop=1) #find radius edge neighbors
            node_ids.extend(neigbhor_node)
            node_ids = set(node_ids)

            subG = G.subgraph(node_ids).copy()
            subG.graph['y'] = peptide.hdx_value
            subG.graph['range'] = (peptide.start, peptide.end)
            subG.graph['chain'] = peptide.chain
            subG.graph['is_complex'] = self.complex_state[complex_state]     

            subG.graph['seq_embedding'] = seq_embedding
            #subG.graph['sequence'] = peptide.sequence
            #subG.graph['pdb'] = apo_identifier

            data = networkx_to_tgG(subG)
            graph_ensemble.append(data)
        label = f'{apo_identifier}'
        return graph_ensemble, label
    
    def _getitem_GVP(self, index):
        '''
        return the graph for the index-th protein complex
        keys format: [database_id, protein, state, match_uni, apo_identifier, chain_identifier, correction, protein_chains]
        '''
        atoms._registry = {}
        res._registry = []
        pep._registry = []

        database_id = self.keys[0][index].strip()
        apo_identifier = self.keys[4][index].strip().split('.')[0].upper()

        protein = self.keys[1][index]
        state = self.keys[2][index]
        match_uni = self.keys[3][index]
        chain_identifier = self.keys[5][index]
        correction = self.keys[6][index]
        protein_chains = self.keys[7][index].split(',')
        complex_state = self.keys[8][index]

        print("processing",database_id, apo_identifier)
        #if os.path.isfile(os.path.join(self.save_dir, f'{apo_identifier}.pt')):
        #    print('skip')
        #    return None

        pdb_fpath = os.path.join(self.pdb_dir, f'{apo_identifier}.pdb')
        target_HDX_fpath = os.path.join(self.hdx_dir, f'{database_id}_revised.xlsx')
        embedding_fpath = os.path.join(self.embedding_dir, f'{apo_identifier}.pt')
        
        if not os.path.isfile(target_HDX_fpath) and self.truncation_window_size is None:
            print("cannot find HDX table: ", target_HDX_fpath)
            return None 
        if not os.path.isfile(embedding_fpath):
            print("Cannot find embedding file:", embedding_fpath)
            return None
        if not os.path.isfile(pdb_fpath):
            print("Cannot find pdb file:", pdb_fpath)
            return None

        # residue graph generation #
        print('creating whole protein graph')
        G = GVP_ResGraph(pdb_fpath, embedding_fpath, protein_chain=protein_chains) # load protein and record backbone coordinations

        residue_coord = np.array([G.nodes[node]['residue_coord'] for node in G.nodes()])
        coords = torch.as_tensor(residue_coord, dtype=torch.float32)
        coords = coords[:,1,:]
        print(coords.shape)
        G = add_edges(G, coord=coords, max_distance=self.max_distance, min_distance=self.min_distance)

        graph_constructor = modified_GVPdataset([], top_k=10)
        G = graph_constructor._featurize_as_graph(G)

        HDX_df = pd.read_excel(target_HDX_fpath, sheet_name='Sheet1')
        ### -------- filter HDX table ------ ####
        HDX_df['state'] = HDX_df['state'].str.replace(' ', '')
        HDX_df['protein'] = HDX_df['protein'].str.replace(' ', '')
        protein = [p.replace(' ', '').strip() for p in protein]
        state = [s.replace(' ', '').strip() for s in state]

        # filter the HDX table
        HDX_df = HDX_df[(HDX_df['state'].isin(state)) & (HDX_df['protein'].isin(protein))]
        HDX_df = HDX_df[HDX_df['sequence'].str.len() < self.max_len]
        HDX_df = HDX_df.sort_values(by=['start', 'end'], ascending=[True, True]) # mixed states will be separated in read_HDX_table func.
        HDX_df = HDX_df.reset_index(drop=True)
        print('after filtering:', HDX_df.shape)
        
        if not read_HDX_table(HDX_df, protein, state, chain_identifier, correction, protein_chains, self.cluster_index): # read HDX table file
            return None
        print('read in pep numbers:', len(pep._registry))

        def find_neigbhors(G, nodes, hop=1):
            neigbhor_node = []
            for node in nodes:
                neigbhors = list(G.neighbors(node))
                neigbhor_node.extend(neigbhors)
            if hop > 1:
                neigbhor_node = find_neigbhors(neigbhor_node, hop-1)
            return neigbhor_node
        
        i2res_id = {(data['chain_id'], data['residue_id']): node for node, data in G.nodes(data=True)}
        graph_ensemble = []
        print('spliting peptide graph')
        for peptide in pep._registry:           
            chain = peptide.chain
            node_ids = [i2res_id[(chain, res_id)] for res_id in peptide.clusters if (chain, res_id) in i2res_id]
            if len(node_ids) == 0:
                continue

            seq_embedding = [G.nodes[node]['x'] for node in node_ids]
            seq_embedding = torch.stack(seq_embedding, dim=0)
            pad_needed = self.max_len - len(node_ids)
            seq_embedding = F.pad(seq_embedding, (0, 0, 0, pad_needed), 'constant', 0)

            neigbhor_node = find_neigbhors(G, node_ids, hop=1)
            node_ids.extend(neigbhor_node)
            node_ids = set(node_ids)

            subG = G.subgraph(node_ids).copy()

            residue_coord = np.array([subG.nodes[node]['residue_coord'] for node in subG.nodes()])
            coords = torch.as_tensor(residue_coord, dtype=torch.float32)
            node_s = torch.stack([subG.nodes[node]['node_s'] for node in subG.nodes()], dim=0)
            node_v = torch.stack([subG.nodes[node]['node_v'] for node in subG.nodes()], dim=0)
            edge_s = torch.stack([subG.edges[(u,v,k)]['edge_s'] for u,v,k in subG.edges(keys=True)], dim=0)
            edge_v = torch.stack([subG.edges[(u,v,k)]['edge_v'] for u,v,k in subG.edges(keys=True)], dim=0)
            mask = torch.isfinite(coords.sum(dim=(1,2)))
            
            Hetero_data = networkx_to_HeteroG(subG)
            Hetero_data['residue'].node_s = node_s
            Hetero_data['residue'].node_v = node_v
            Hetero_data['residue'].edge_s = edge_s
            Hetero_data['residue'].edge_v = edge_v
            Hetero_data['residue'].mask = mask

            '''data = torch_geometric.data.Data(y = torch.as_tensor(peptide.hdx_value), range = torch.as_tensor([peptide.start, peptide.end]), chain = torch.as_tensor(peptide.chain), 
                                             seq_embedding = seq_embedding, is_complex = torch.as_tensor(self.complex_state[complex_state]),
                                            node_s=node_s, node_v=node_v,
                                            edge_s=edge_s, edge_v=edge_v,
                                            edge_index=edge_index_mapped, mask=mask)'''
            
            graph_ensemble.append(Hetero_data)

        label = f'{apo_identifier}'
        return graph_ensemble, label