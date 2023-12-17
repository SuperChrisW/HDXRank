#################################################################################################
# This script is used to calculate the protein rigidity and flexibility based on the ANSSURR method
# The ANSSURR method is based on the FIRST method, which is a rigidity index for each residue
# WANG LIYAO, 2023-12-17
#################################################################################################

import os
import shutil
import numpy as np
from collections import Counter
from Bio.PDB import PDBParser, is_aa
from Bio.PDB.DSSP import make_dssp_dict, DSSP
import warnings

#from biotite.structure import AtomArray, AtomArrayStack, sasa
#import biotite.structure.io.pdb as pdb

warnings.filterwarnings("ignore")

################################# ANSSURR rigidity process block #################################
class res(object):
    _registry = []
    def __init__(self, i, name, atom):
        self._registry.append(self)
        self.clusters = []
        self.energy = 0
        self.i = i
        self.name = name
        self.CA = atom

def rescale_FIRST(FIRST, smooth_ = True):
    K = 315777.09
    T = 298.15
    Hartree = 0.00038
    result = np.exp((4.2 * FIRST * Hartree) * K / T)
    result[np.isnan(FIRST)] = np.nan
    #transverse result matrix
    if smooth_:    
        result = result.T
        for i in range(len(result)):
            result[i] = smooth([res.i for res in res._registry], result[i])
        return result.T
    else:
        return result

def smooth(resi,data):
    data_smoothed = []
    segs = get_segments(resi)
    for s in segs:
        if len(s) > 2:
            data_smoothed_temp = movingaverage(data[s[0]:s[1]+1],3)
            N = (data[s[0]] + data[s[0]+1]) / 2.0
            data_smoothed_temp[0] = N
            C = (data[s[1]] + data[s[1]-1]) / 2.0
            data_smoothed_temp[-1] = C
            data_smoothed.extend(data_smoothed_temp)
        elif len(s) == 2:
            data_smoothed_temp = movingaverage(data[s[0]:s[1]+1],2)
            N = (data[s[0]] + data[s[0]+1]) / 2.0
            data_smoothed_temp[0] = N
            C = (data[s[1]] + data[s[1]-1]) / 2.0
            data_smoothed_temp[-1] = C
            data_smoothed.extend(data_smoothed_temp)
        else:
            data_smoothed.append(data[s[0]])
    return data_smoothed

def movingaverage(interval, window_size):
    window= np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')

def get_segments(resi):
    prev_res = min(resi) if resi else None
    segments = list()
    for number in enumerate(resi):
        if number[1] != prev_res+1:
            segments.append([number[0]])
        elif len(segments[-1]) > 1:
            segments[-1][-1] = number[0]
        else:
            segments[-1].append(number[0])
        prev_res = number[1]
    return segments

def index_consecutive_letters(letter_list):
    indices = {}
    current_letter = letter_list[0]
    ss_count = 0

    for i, letter in enumerate(letter_list):
        if letter == current_letter:
            if ss_count not in indices.keys():
                indices[ss_count] = [i]
            else:
                indices[ss_count].append(i)
        else:
            ss_count += 1
            indices[ss_count] = [i]
            current_letter = letter
    return indices

###############################################################################################



################################# protein property output block #################################
### output ASSURR residue-wise flexibility ###
def ASSURR_flexibility(proflex_path, dssp_path, save_path):
## FIXME: the file path is not correct, need to be fixed
    
    if not os.path.exists(proflex_path) and not os.path.exists(dssp_path):
        print('cannot find the folder')
        return False
    for file in os.listdir(dssp_path):
        PDB = file.split('.')[0]
        if PDB == '': continue
        if os.path.exists(f'{save_path}/rigidity_{PDB}.txt'): continue

        ss_structure = []        
        dssp_RSA_path = os.path.join(dssp_path, f'{PDB}.dssp.txt')
        with open(dssp_RSA_path, 'r') as f:
            data = f.readlines()
            for line in data:
                ss_value = line.split()[2]
                ss_structure.append(ss_value)
        ss_structure = np.array(ss_structure)

        energy_list_m = []
        #try:     
        print(f'processing {PDB}')   
        linkref, ph_bonds, torsions, CA_mut, natom, _, _ = read_edge(f'{PDB}_proflexdataset', f'{proflex_path}', 0, central_atom= 'CA')
        rigid_clst_m, clst_ensemble_m = read_decomp_list(f'decomp_{PDB}', f'{proflex_path}', natom, 0, ss_list=ss_structure)
        for residue in res._registry:
            energy_list_m.append(residue.energy)
        
        energy_list_m = np.reshape(energy_list_m, (len(energy_list_m), 1))

        rescale_rigid_m = rescale_FIRST(energy_list_m, smooth_ = True)
        np.savetxt(f'{save_path}/rigidity_{PDB}.txt', rescale_rigid_m)
        #except Exception as e:
        #    print(f'cannot process {PDB} as {e}')
        #    continue
    return True

### output relative solvent accessible surface area (SASA) by biopython.PDB.DSSP ###
### The solvent accessible surface area (SASA) was calculated using BioPython and default Sander and Rost values (Cock et al., 2009; Rost and Sander, 1994). 

def DSSP_calculation(pdb_folder_path, save_dir):
    if not os.path.exists(pdb_folder_path):
        print('cannot find the pdb folder: ', pdb_folder_path)
        return False
    for file in os.listdir(pdb_folder_path):
        if not file.endswith('.pdb'):
            continue
        file = file.replace('.pdb', '')
        PDB_path = os.path.join(f'{pdb_folder_path}', f'{file}.pdb')

        parser = PDBParser()
        structure = parser.get_structure(file, PDB_path)

        model = structure[0]  # Assuming there's only one model in the structure
        dssp = DSSP(model, PDB_path)

        # calculate the number of AA in the protein 
        aa_list = []
        for chain in structure.get_chains():
            AA_len = len([_ for _ in chain.get_residues() if is_aa(_)])
            
            for id in chain.get_residues():
                if is_aa(id):
                    aa_list.append(id)

        out_path = os.path.join(save_dir, f'{file}.dssp.txt')
        with open(out_path,"w") as out:
            for i in range(0,len(aa_list)):
                a_key = list(dssp.keys())[i]
                for item in list(dssp[a_key]):
                    out.write(str(item)+"\t")
                out.write("\n")
        print(f'finished {file}')


### output relative solvent accessible surface area (SASA) by biotite ###
### The solvent accessible surface area (SASA) was calculated in atom-wise.

max_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0
} # TIEN et al., 2013 (theoretical maximum solvent accessible surface area)

def biotite_SASA(file_list, root_dir, fname): ## FIXME: the file path is not correct, need to be fixed
    for file in file_list:
        PDB = file.split('.')[0]
        print(PDB)
        if PDB == '': continue
        #if PDB != '1OPA_wiH': continue
        PDB_folder = os.path.join(f'{root_dir}/HDX_Dynamic', f'{PDB}')

        ensemble = [f'{PDB}{suffix}' for suffix in fname]
        biotite_RSA = []
        for ens in ensemble:
            RSA_N = []
            src_path = os.path.join(PDB_folder, ens)
            if os.path.exists(src_path):
                structure = pdb.PDBFile.read(src_path)
                model = structure.get_structure(model=1)
                N_mask = model.atom_name == "N"
                H_mask = model.atom_name == "H"

                SASA_atom = sasa(model, atom_filter = N_mask, point_number= 100, vdw_radii="ProtOr")
                SASA_H = sasa(model, atom_filter = H_mask, point_number= 100, vdw_radii="Single")

                for i in range(len(SASA_atom)):
                    if N_mask[i] and H_mask[i+1]: 
                        #print(N_mask[i], H_mask[i+1], i)

                        #print(SASA_atom[i], SASA_H[i+1])
                        SASA_atom[i] += SASA_H[i+1]

                not_nan_mask = ~np.isnan(SASA_atom)
                SASA_N = SASA_atom[not_nan_mask]
                res_name = model.res_name[N_mask]

                for asa, res in zip(SASA_N, res_name):
                    max_asa = max_ASA.get(res)
                    if max_asa is not None:
                        #rsa = asa / max_asa
                        rsa = asa
                        RSA_N.append(rsa)
                        #print('residue: ', res, asa, max_asa, 'rsa: ', rsa)
                    else:
                        print(f"Warning: Max ASA value not found for residue {res}")
                        RSA_N.append(np.nan)
                biotite_RSA.append(RSA_N)

        biotite_RSA = np.stack(biotite_RSA, axis=0)
        avg_biotite_RSA = np.mean(biotite_RSA, axis=0)

        with open(f'{PDB_folder}/biotite_SASA_test{PDB}.txt', 'w') as f:
            for i in range(len(avg_biotite_RSA)):
                f.write(f'{i} {avg_biotite_RSA[i]}\n')
        print(f'finished {PDB}')
        break
###############################################################################################


################################ proflex output file process block ############################
    
def read_edge(key, direct, H_bond_cutoff = -1.0, central_atom = 'CA'): # read the proflex dataset file and return the bonding information as edge_list
    path = os.path.join(direct,key)
    nhb, nph = -1, -1 # index of hbonds, hydrophobic interactions

    clst = {} # node dict = {node:[#atoms]}
    rev_clst = {} # reverse node dict = {atom: #node}
    linkref = {} #create linkref list to represent the neighbor atoms/ adj matrix for atoms, including non-covalent bonding
    hbonds = {} # store the informaiton 
    H_energy = {} # store the bond energy
    bond_type = {} # store the non-covalent bond type

    cf_list = [] # central-force bond list
    ph_bonds = [] # hydrophobic interaction , ph_bonds[index] = [[real atom, 3 psedo-atoms index]]
    torsions = [] # record the locked central force pairs
    CA_list = {} # record the CA index
    res._registry = [] # reset the residue registry
    natom, nres= 0, 0 
    descritpion = {} # store the description of H bond

#    print(path)
    if not os.path.isfile(path):
        print("cannot find the file", key)
        return None

    with open(path, 'r') as f:
        data = f.read().strip().split('\n') 
        for line in data:
            ############## clustering atoms into residue group ###############
            if line[:4] == 'ATOM':
                resn = line[17:20].replace(' ','') # resn is residue name, remove spaces
                natom += 1 
                res_id = int(line[22:26])
                if line[13:15].strip() == central_atom: ## should be CA
                    CA_list[res_id] = int(line[6:11])-1 # CA index starts from 0 
                    residues = res(res_id, resn, natom)

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

            for bond in ph_bonds:
                so = bond[0]
                sf = bond[-1]
                
    #    print("total residue: ", len(clst.keys()))
    #    print("total atom: ", len(rev_clst.keys()))
#    print(nph)
    return linkref, ph_bonds, torsions, CA_list, natom, (hbonds, H_energy, bond_type), descritpion

def read_decomp_list(filename, directory, natom, cutoff = -1.0, ss_list = None):
    compressed_data = []  # atom to rigid cluster mapping 
    current_array = []  # Temporary list for one read operation
    decomp_cluster = {} # rigid cluster to atom mapping
    rigid_clusters = [] # rigid cluster list (contain atoms >= 15)
    H_energy = []

    if len(ss_list) > 0:
        ss_clst = index_consecutive_letters(ss_list)
        ss_res = [key for key, value in ss_clst.items() for _ in value]

    # Open the file for reading
    with open(f"{directory}/{filename}", 'r') as file:
        skip = True  # Initialize skip to True
        for line in file:
            line = line.strip()
            if line.startswith(("A", "B", "HEADER")):
                parts = line.split()
                #print(parts)
                if parts[0] == 'A:':
                    energy = float(parts[2])
                    if energy <= cutoff:
                        skip = False
                        H_energy.append(energy)
                    else:
                        skip = True
                continue
            elif line.startswith("END"):
                #if compressed_data: # only take the rigid cluster at transition state 
                #    break
                if not skip and current_array:
                    rc = []
                    current_array = np.array(current_array[:natom])
                    compressed_data.append(current_array)

                    cluster_count = Counter(current_array)

                    for c in cluster_count:
                        if int(cluster_count[c]) >= 15: # c is cluster index, if it contains more than 15 atoms (about 2 residues), append to rigid cluster list
                            rc.append(c)
                    rigid_clusters.append(rc)
                    current_array = []
            elif not skip: 
                try:
                    line = line.replace('****', '9999')
                    numbers = [int(num) for num in line.split(':') if num.strip()]
                    current_array.extend(numbers)
                except ValueError:
                    print(line)
                    pass

    #clst_ensemble = np.column_stack(compressed_data)
    clst_ensemble = np.array(compressed_data)

    for d in clst_ensemble:
        for r in res._registry: # go through each residue, it stores CA atom info for each residue
            r.clusters.append(d[r.CA-1]) # r.CA-1 is the atom number, d[r.CA-1] is the cluster index for that atom
            #print(r.name, r.i, r.CA, r.clusters[-1])
        #print(r.clusters) # r is the last residue in the registry

    for x in range(len(clst_ensemble)):
        for r in res._registry:
            if r.clusters[x] in rigid_clusters[x]: # if this res is belong to a rigid cluster, add energy
                
                res_id = r.i
                ss_type = ss_list[res_id-1]
                ss_clst_id = ss_res[res_id-1]
                ss_clst_num = len(ss_clst[ss_clst_id])
                #print(ss_type, ss_clst_num)
                
                if ss_type in ['H', 'G', 'I', 'E']:
                    if (ss_type == 'H' and ss_clst_num > 4) or (ss_type == 'G' and ss_clst_num > 3) or \
                        (ss_type == 'I' and ss_clst_num > 5) or (ss_type == 'E' and ss_clst_num > 3):
                        r.energy = H_energy[x]
                else:
                    r.energy = H_energy[x]

                #r.energy = H_energy[x] # r.energy records the last energy value for that residue in rigid cluster

    return rigid_clusters, clst_ensemble

############################################################################################################