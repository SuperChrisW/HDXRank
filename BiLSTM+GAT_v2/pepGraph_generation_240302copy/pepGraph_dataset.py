import numpy as np
import pandas as pd

from scipy.spatial import distance_matrix
import torch
from torch.utils.data import Dataset
from torch.utils.data.sampler import Sampler

import os
from Bio.PDB.Polypeptide import protein_letters_3to1
from pepGraph_utlis import load_embedding

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
    def __init__(self, i, name, atom, chain_id, position):
        self._registry.append(self)
        self.clusters = []
        self.energy = 0
        self.i = i
        self.name = name
        self.CA = atom
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
    protein_embedding = {}
    def __init__(self, i, pep_sequence, chain, start_pos, end_pos, hdx_value, log_t):
        self._registry.append(self)
        self.i = i
        self.sequence = pep_sequence
        self.chain = chain
        self.start = start_pos
        self.end = end_pos
        self.clusters = []
        self.position = None
        self.embedding = None
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
        return adj
    
    @classmethod
    def get_pep_feature(cls, embedding_fpath): # obtain the feature matrix for peptide
        ### Note: embedding_fpath is a list of embedding file paths for each chain, {apo_identifier}_{chain}.embedding.txt
        for fpath in embedding_fpath:
            current_chain = fpath.split('.')[0].split('_')[-1]
            embedding_dict, embedding_matrix, _ = load_embedding(fpath)

            pep.protein_embedding[current_chain] = embedding_matrix
            for peptide in cls._registry:
                start = peptide.start
                end = peptide.end
                chain = peptide.chain
                if chain == current_chain:
                    mtx = torch.zeros((end-start+1, embedding_matrix.shape[1]), dtype = torch.float32)
                    for index, res_id in enumerate(range(start-1, end)):
                        if res_id in embedding_dict.keys():
                            mtx[index] = embedding_dict[res_id].to(dtype = torch.float32)

                    rigidity_std = np.std(mtx[:, 1].view(-1).numpy())
                    seq_length = mtx.shape[0] / 30
                    extra_vec = (peptide.log_t, rigidity_std, seq_length)
                    extra_col = torch.Tensor(extra_vec)
                    mtx_mean = torch.mean(mtx, dim = 0)
                    peptide.embedding = torch.cat((mtx_mean, extra_col), dim = 0)
                    
                    if np.isnan(peptide.embedding).any():
                        print(peptide.embedding)                        
                        peptide.embedding[torch.isnan(peptide.embedding)] = 0
                        print(f"Embedding matrix contains NaN values: {fpath}")
                        break


    @classmethod
    def set_pep_position(cls): # obtain the average position of peptide
        for peptide in cls._registry:
            pep_position = np.array([0,0,0], dtype = np.float32)
            for res_index in peptide.clusters:
                residue = res.get_res(res_index, peptide.chain)
                pep_position += residue.position

            num_residues = len(peptide.clusters)
            peptide.position = pep_position/num_residues

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
    max_len = 30
    for residue in res._registry:
        if residue.chain_id not in res_chainsplit.keys():
            res_chainsplit[residue.chain_id] = []
        res_chainsplit[residue.chain_id].append(residue)

    ## seperately process states
    for index, (protein, state, chain) in enumerate(zip(proteins, states, chains)):
        temp_HDX_df = HDX_df[(HDX_df['state']==state) & (HDX_df['protein']==protein)]
        temp_HDX_df = temp_HDX_df[temp_HDX_df['sequence'].str.len() < max_len]
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
        pep.set_pep_position()
        return True

class pepGraph(Dataset):
    def __init__(self, keys, root_dir, max_len=30):
        self.keys = keys
        self.root_dir = root_dir
        self.embedding_dir = os.path.join(root_dir, 'embedding_files')
        self.proflex_dir = os.path.join(root_dir, 'proflex_files')
        self.hdx_dir = os.path.join(root_dir, 'HDX_files')
        self.max_len = max_len

    def __len__(self):
        return len(self.keys)

    def __getitem__(self, index):
        atoms._registry = {}
        res._registry = []
        pep._registry = []
        pep.protein_embedding = {}

        database_id = self.keys[0][index].strip()
        apo_identifier = self.keys[4][index].strip().split('.')[0]

        protein = self.keys[1][index]
        state = self.keys[2][index]
        match_uni = self.keys[3][index]
        chain_identifier = self.keys[5][index]
        correction = self.keys[6][index]
        print("processing",database_id, apo_identifier)
            
        proflex_fpath = os.path.join(self.proflex_dir, f'{apo_identifier}_clean_Hplus_proflexdataset')
        target_HDX_fpath = os.path.join(self.hdx_dir, f'{database_id}.xlsx') 
        embedding_fpath = [os.path.join(self.embedding_dir, f'{apo_identifier}_{chain}.embedding.txt') for chain in chain_identifier]
        
        if not os.path.isfile(proflex_fpath):
            print("cannot find proflex file: ", proflex_fpath)
            return None 
        elif not os.path.isfile(target_HDX_fpath):
            print("cannot find HDX table: ", target_HDX_fpath)
            return None 
        for path in embedding_fpath:
            if not os.path.isfile(path):
                print("Cannot find embedding file:", path)
                return None

        HDX_df = pd.read_excel(target_HDX_fpath, sheet_name='Sheet1')
        ### -------- filter HDX table ------ ####
        HDX_df['state'] = HDX_df['state'].str.replace(' ', '')
        HDX_df['protein'] = HDX_df['protein'].str.replace(' ', '')

        protein = [p.replace(' ', '').strip() for p in protein]
        state = [s.replace(' ', '').strip() for s in state]

        #print(HDX_df['protein'].unique(), HDX_df['state'].unique())
        #print(protein, state)
        HDX_df = HDX_df[HDX_df['state'].isin(state)]
        HDX_df = HDX_df[HDX_df['protein'].isin(protein)]

        HDX_df = HDX_df.drop_duplicates(subset=['sequence', 'state', 'protein'])
        HDX_df = HDX_df[HDX_df['sequence'].str.len() < self.max_len]

        HDX_df = HDX_df.sort_values(by='start', ascending=True) # mixed states will be separated in read_HDX_table func.
        HDX_df = HDX_df.reset_index(drop=True)
        print('after filtering:', HDX_df.shape)

        ### -------- graph generation ------ #### 
        linkref, ph_bonds, torsions, CA_list, natom, (hbonds, H_energy, bond_type), descritpion = \
            read_edge(apo_identifier, proflex_fpath, H_bond_cutoff= -1.0, central_atom= 'CA') # read proflex dataset file
        print("number of atoms:", len(atoms._registry.keys()))
        print('max index of atom:', max(list(linkref.keys())))
        print("number of residues:", len(res._registry))

        if not read_HDX_table(HDX_df, protein, state, chain_identifier, correction): # read HDX table file
            return None
        print('read in pep numbers:', len(pep._registry))

        '''
        for atom in atoms._registry.values():
            print(atom.i, atom.name, atom.res_cluster)
        for residue in res._registry:
            print(residue.i, residue.name, residue.clusters)
        '''

        pep_adj = pep.get_pep_adj(linkref, ph_bonds) # obtain the adjacency matrix for peptide

        I = np.eye(len(pep._registry))
        pep_adj = pep_adj + I # node self-connection
        pep.get_pep_feature(embedding_fpath)

        label = f'{apo_identifier}'
        return label, pep._registry, pep_adj
