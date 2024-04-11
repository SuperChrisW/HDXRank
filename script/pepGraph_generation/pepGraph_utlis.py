import os
import re
import math
import pandas as pd
import numpy as np
import torch
import torch.nn.functional as F
from pdb2sql import interface
from tool.BioWrappers import *
from dataclasses import dataclass, field
from Bio.PDB.Polypeptide import one_to_three as protein_letters_1to3, three_to_one as protein_letters_3to1
from Bio.PDB import Selection
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel

import warnings
warnings.filterwarnings("ignore")

class ChemData():
    def __init__(self):
        self.NAATOKENS = 20+1+10+10+1 # 20 AAs, 1 UNK res, 8 NAs+2UN NAs, 10 atoms +1 UNK atom
        self.UNKINDEX = 20  # residue unknown

        #bond types
        self.num2btype = [0,1,2,3,4,5,6,7] # UNK, SINGLE, DOUBLE, TRIPLE, AROMATIC, 
                                            # PEPTIDE/NA BACKBONE, PROTEIN-LIGAND (PEPTIDE), OTHER

        self.NATYPES = ['DA','DC','DG','DT', 'DX', 'RA','RC','RG','RU', 'RX']
        self.STDAAS = ['ALA','ARG','ASN','ASP','CYS',
            'GLN','GLU','GLY','HIS','ILE',
            'LEU','LYS','MET','PHE','PRO',
            'SER','THR','TRP','TYR','VAL',]
        
        self.num2aa=[
            'ALA','ARG','ASN','ASP','CYS',
            'GLN','GLU','GLY','HIS','ILE',
            'LEU','LYS','MET','PHE','PRO',
            'SER','THR','TRP','TYR','VAL',
            'UNK',
            'DA','DC','DG','DT', 'DX',
            'RA','RC','RG','RU', 'RX', 
            'Br', 'F', 'Cl','I',
            'C', 'N', 'O', 'S', 'P', 'Si',
            'ATM'
        ]
        self.num2aa = [item.upper() for item in self.num2aa]
        self.aa2num= {x:i for i,x in enumerate(self.num2aa)}
        self.aa2num['MEN'] = 20

        # Mapping 3 letter AA to 1 letter AA (e.g. ALA to A)
        self.one_letter = ["A", "R", "N", "D", "C", \
                            "Q", "E", "G", "H", "I", \
                            "L", "K", "M", "F", "P", \
                            "S", "T", "W", "Y", "V", "?", "-"]

        self.n_non_protein = len(self.num2aa) - len(self.one_letter)

        self.aa_321 = {a:b for a,b in zip(self.num2aa,self.one_letter+['a']*self.n_non_protein)}

        self.frame_priority2atom = [
            "F",  "Cl", "Br", "I",  "O",  "S",  "Se", "Te", "N",  "P",  "As", "Sb", 
            "C",  "Si", "Sn", "Pb", "B",  "Al", "Zn", "Hg", "Cu", "Au", "Ni", "Pd", 
            "Pt", "Co", "Rh", "Ir", "Pr", "Fe", "Ru", "Os", "Mn", "Re", "Cr", "Mo", 
            "W",  "V",  "U",  "Tb", "Y",  "Be", "Mg", "Ca", "Li", "K",  "ATM"]

        # these atomic numbers are incorrect, but keeping for fold&dock3 and correcting it 
        # in util.writepdb() during output.
        self.atom_num= [
            9,    17,   35,   53,   8,    16,   34,   52,   7,    15,   33,   51, 
            6,    14,   32,   50,   82,   5,    13,   30,   80,   29,   79,   28, 
            46,   78,   27,   45,   77,   26,   44,   76,   25,   75,   24,   42, 
            23,   74,   92,   65,   39,   4,    12,   20,   3,    19,   0] # in same order as frame priority

        self.atom2frame_priority = {x:i for i,x in enumerate(self.frame_priority2atom)}
        self.atomnum2atomtype = dict(zip(self.atom_num, self.frame_priority2atom))
        self.atomtype2atomnum = {v:k for k,v in self.atomnum2atomtype.items()}
                
        self.residue_charge = {'CYS': -0.64, 'HIS': -0.29, 'ASN': -1.22, 'GLN': -1.22, 'SER': -0.80, 'THR': -0.80, 'TYR': -0.80,
                                'TRP': -0.79, 'ALA': -0.37, 'PHE': -0.37, 'GLY': -0.37, 'ILE': -0.37, 'VAL': -0.37, 'MET': -0.37,
                                'PRO': 0.0, 'LEU': -0.37, 'GLU': -1.37, 'ASP': -1.37, 'LYS': -0.36, 'ARG': -1.65}

        self.residue_polarity = {'CYS': 'polar', 'HIS': 'polar', 'ASN': 'polar', 'GLN': 'polar', 'SER': 'polar', 'THR': 'polar', 'TYR': 'polar', 'TRP': 'polar',
                                    'ALA': 'apolar', 'PHE': 'apolar', 'GLY': 'apolar', 'ILE': 'apolar', 'VAL': 'apolar', 'MET': 'apolar', 'PRO': 'apolar', 'LEU': 'apolar',
                                    'GLU': 'neg_charged', 'ASP': 'neg_charged', 'LYS': 'neg_charged', 'ARG': 'pos_charged'}

        self.polarity_encoding = {'apolar': 0, 'polar': 1, 'neg_charged': 2, 'pos_charged': 3}

        self.ss_list = ['H', 'B', 'E', 'G', 'I', 'T', 'S', 'P', '-']

        self.AA_array = {
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
chemdata = ChemData()

@dataclass
class RawInputData:
    msa: torch.Tensor = field(default_factory=torch.Tensor)
    res_HDMD: torch.Tensor = field(default_factory=torch.Tensor)
    res_polarity: torch.Tensor = field(default_factory=torch.Tensor)
    res_charge: torch.Tensor = field(default_factory=torch.Tensor)
    SASA: torch.Tensor = field(default_factory=torch.Tensor)
    rigidity: torch.Tensor = field(default_factory=torch.Tensor)
    hse: torch.Tensor = field(default_factory=torch.Tensor)
    seq_data: dict = field(default_factory=dict)  # [id]: {'token_type', 'coord'}
    type_label: str = ''

    def construct_embedding(self):
        # construct embedding
        embedding = torch.cat((self.msa, self.res_HDMD, self.res_polarity, self.res_charge, self.SASA, self.hse), dim=1)
        # 240410 delete self.rigidity from the embedding
        return embedding
    
    def merge(self, data):
        self.seq_data = {**self.seq_data, **data.seq_data}
        return self

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

    return embedding_dict, embedding_matrix, embedding_seq

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

    start_pos = [x + correction_value for x in df['start'].astype(int).tolist()]
    end_pos = [x + correction_value for x in df['end'].astype(int).tolist()]
    seq = df['sequence'].tolist()
    log_t = df['log_t'].to_numpy()

    # read embedding file
    embed_dict, embed_array, embed_seq = load_embedding(vector_file, feat_rescale=False)
    embed_seq = [protein_letters_3to1(res) for res in embed_seq]

    row_size = len(start_pos)
    nfactor = embed_array.shape[-1] + 3
    input_array = np.zeros((row_size,max_len,nfactor), dtype = np.float32)

    # filter out the sequences that do not match the embedding file
    pop_idx = []

    for i, (start, end, sequence) in enumerate(zip(start_pos, end_pos, seq)):
        embed_sequence = ''.join(embed_seq[start - 1:end])
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

# feature extraction
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
    encode_index = [chemdata.polarity_encoding[chemdata.residue_polarity[protein_letters_1to3(res)]] for res in seq]
    polarity_mtx = np.zeros((len(seq), 4))
    for i, idx in enumerate(encode_index):
        polarity_mtx[i, idx] = 1
    return polarity_mtx

def parse_hhm(hhm_file):
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

# loading protine, NA, molecule ligand
def load_protein(hhm_file, rigidity_file, pdb_file, chain_id):
    '''
    Generate embedding file from the pre-computed HMM file and rigidity file, PDB strucutre
    processing the protein chain featurization
    '''
    ### convert the list in dictionary to array

    model = get_bio_model(pdb_file)
    residue_data = {}
    res_seq = []
    residue_list = Selection.unfold_entities(model, 'R')
    for res in residue_list:
        if res.get_parent().get_id() != chain_id:
            continue
        #if res.get_id()[0] != ' ': # exclude HETATM
        #    continue
        res_id = res.get_id()[1]
        res_name = res.get_resname()
        res_seq.append(protein_letters_3to1(res_name) if res_name in chemdata.STDAAS else 'X')

        ca_atom = res['CA'] if 'CA' in res else None
        if ca_atom is None: raise ValueError(f'CA atom not found in the residue at {res_id}')
        ca_coord = ca_atom.get_coord()
        residue_data[res_id] = {
            'token_type': res_name,
            'coord': ca_coord,
        }

    rigidity = torch.tensor([])
    #rigidity_df = pd.read_table(rigidity_file, sep=' ', header=None)
    #rigidity_df.columns = ['residue_index', 'residue_name', 'chain_index', 'rigidity']
    #rigidity_df = rigidity_df[rigidity_df['chain_index'] == chain_id]

    max_len = len(res_seq)
    #rigidity_df = rigidity_df.iloc[:max_len]
 
    # sequence label
    #res_id = rigidity_df['residue_index'].to_numpy().reshape(-1, 1)
    #res_name = rigidity_df['residue_name'].to_numpy().reshape(-1, 1)

    # sequence-based features
    res_charge = np.array([chemdata.residue_charge[protein_letters_1to3(res)] for res in res_seq]).reshape(-1, 1)
    res_polarity = get_seq_polarity(res_seq)
    HDMD = np.array([chemdata.AA_array[res] for res in res_seq]).reshape(-1, 5)

    # structure-based features
    hse_dict = get_hse(model, chain_id)
    SASA = biotite_SASA(pdb_file, chain_id)[:max_len]
    #rigidity = rigidity_df['rigidity'].to_numpy().reshape(-1, 1)

    # MSA-based features
    hhm_mtx, hhm_seq = parse_hhm(hhm_file) # hhm file is chain-wise

    # check the sequence match among feat.
    res_seq = ''.join(res_seq)
    if hhm_seq == res_seq:
        pass
    elif hhm_seq[:-1] == res_seq:
        hhm_mtx = hhm_mtx[:-1]
    else:
        print(hhm_seq)
        print(res_seq)
        raise ValueError('Sequence mismatch between HMM and DSSP')

    corrected_hse_mtx = np.zeros((max_len, 3))
    for i, res_j in enumerate(residue_data.keys()):
        res_j = str(res_j)
        if (chain_id, res_j) in hse_dict.keys():
            corrected_hse_mtx[i, :] = list(hse_dict[(chain_id, res_j)])

    print('protein length:', len(residue_data.keys()))
    #print('rigidity length:', rigidity.shape)
    print('dssp length:', SASA.shape)
    print('hse length:', corrected_hse_mtx.shape)
    print('HDMD length:', HDMD.shape)
    print('res_charge length:', res_charge.shape)
    print('res_polarity length:', res_polarity.shape)
    print('hhm length:', hhm_mtx.shape)

    return RawInputData(
        msa = torch.tensor(hhm_mtx, dtype = torch.float32),
        res_HDMD = torch.tensor(HDMD, dtype = torch.float32),
        res_polarity = torch.tensor(res_polarity, dtype = torch.float32),
        res_charge = torch.tensor(res_charge, dtype = torch.float32),
        SASA = torch.tensor(SASA, dtype = torch.float32).reshape(-1, 1),
        rigidity = torch.tensor(rigidity, dtype = torch.float32),
        hse = torch.tensor(corrected_hse_mtx, dtype = torch.float32),
        seq_data = residue_data,
        type_label = 'protein'
    )

def load_sm(pdb_file, preset_chain_id):
    atom_data = {}
    structure = get_bio_model(pdb_file)
    with open(pdb_file, 'r') as f:
        data = f.read().strip().split('\n')
        for line in data:
            if line[:6] == 'HETATM':
                chain_id = line[21]
                if chain_id != preset_chain_id:
                    continue
                LG_name = line[17:20].replace(' ','') # resn is residue name, remove spaces
                atom_id = int(line[6:11].strip())
                atom_type = line[12:16].strip() # only allows for one character symbol for atom type
                element_symbol_regex = ''.join(re.findall('[A-Za-z]', atom_type)).upper()
                token_type = chemdata.aa2num[element_symbol_regex] if element_symbol_regex in chemdata.aa2num else chemdata.aa2num['ATM']
                res_id = line[22:26]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atom_data[atom_id] = {
                    'token_type': token_type,
                    'coord': [x, y, z]
                }

        if len(atom_data.keys()) == 0:
            print('No small molecule found in the PDB file')
            return None
    return RawInputData(
        msa = None,
        res_HDMD = None,
        res_polarity = None,
        res_charge = None,
        SASA = None,
        rigidity = None,
        hse = None,
        seq_data = atom_data,
        type_label = 'ligand'
    )

def load_nucleic_acid(pdb_file, chain_id):
    na_data = {}
    model = get_bio_model(pdb_file)

    for chain in model.get_chains():
        if chain.get_id() != chain_id:
            continue
        for res in chain:
            na_id = res.get_id()[1]
            na_name = res.get_resname().strip()

            if na_name not in chemdata.NATYPES:
                continue

            atom_coord = [atom.get_coord() for atom in res.get_atoms()]
            na_coord = np.mean(atom_coord, axis=0)
            na_data[na_id] = {
                'token_type': chemdata.aa2num[na_name],
                'coord': na_coord,
            }

    return RawInputData(
        msa = None,
        res_HDMD = None,
        res_polarity = None,
        res_charge = None,
        SASA = None,
        rigidity = None,
        hse = None,
        seq_data = na_data,
        type_label = 'NA'
    )

'''
def extract_sm_to_sdf(pdb_file, chain_id, save_dir):
    # extract small molecule from PDB file
    # convert to sdf file
    atoms_info = []
    with open(pdb_file, 'r') as f:
        data = f.read().strip().split('\n')
        for line in data:
            if line[:6] == 'HETATM':
                chain_id = line[21]
                if chain_id != chain_id:
                    continue
                LG_name = line[17:20].replace(' ','') # resn is residue name, remove spaces
                atom_id = int(line[6:11].strip())
                atom_type = line[12:16].strip() # only allows for one character symbol for atom type
                res_id = line[22:26]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms_info.append((atom_type, atom_id,  x, y, z))
        if len(atoms_info) == 0:
            print('No small molecule found in the PDB file')
            return None
            
    mol = Chem.RWMol()
    conf = Chem.Conformer(len(atoms_info))

    for element, _, x, y, z in atoms_info:
        atom = Chem.Atom(element)
        idx = mol.AddAtom(atom)
        conf.SetAtomPosition(idx, (x, y, z))
    mol.AddConformer(conf)

    fname = os.path.basename(pdb_file).split('.')[0]
    fpath = os.path.join(save_dir, f'{fname}_ligand.sdf')
    w = Chem.SDWriter(fpath)
    w.write(mol)
    w.close()
    return fpath, atoms_info

def clean_sdffile(filename):
    # lowercase the 2nd letter of the element name (e.g. FE->Fe) so openbabel can parse it correctly
    lines2 = []
    with open(filename) as f:
        lines = f.readlines()
        num_atoms = int(lines[3][:3])
        for i in range(len(lines)):
            if i>=4 and i<4+num_atoms:
                lines2.append(lines[i][:32]+lines[i][32].lower()+lines[i][33:])
            else:
                lines2.append(lines[i])
    molstring = ''.join(lines2)

    return molstring

def load_sm(pdb_file, chain_id):
    fpath, atoms_info = extract_sm_to_sdf(pdb_file, chain_id, './')    
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat('sdf')
    obmol = openbabel.OBMol()
    molstring = clean_sdffile(fpath)

    if obConversion.ReadString(obmol,molstring):
        print("Molecule loaded successfully.")
    else:
        print("Failed to load the molecule.")
    
    atom_data = {}
    for i in range(1, obmol.NumAtoms()+1):
        atom = obmol.GetAtom(i)
        atom_type = ChemData().atomnum2atomtype.get(atom.GetAtomicNum(), 'ATM')
        atom_data[i] = {
            'token_type': ChemData().aa2num[atom_type],
            'coord': [atom.x(), atom.y(), atom.z()]
        }
    os.remove(fpath) 

    atomtypes = [ChemData().atomnum2atomtype.get(obmol.GetAtom(i).GetAtomicNum(), 'ATM') 
                for i in range(1, obmol.NumAtoms()+1)]
    msa = torch.tensor([ChemData().aa2num[x] for x in atomtypes])
    atom_coords = torch.tensor([[obmol.GetAtom(i).x(),obmol.GetAtom(i).y(), obmol.GetAtom(i).z()]
                                for i in range(1, obmol.NumAtoms()+1)]).unsqueeze(0) # (1, natoms, 3)

    return RawInputData(
        msa = None,
        res_HDMD = None,
        res_polarity = None,
        res_charge = None,
        SASA = None,
        rigidity = None,
        hse = None,
        seq_data = atoms_info,
        type_label = 'ligand'
    )
'''

# graph modification
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