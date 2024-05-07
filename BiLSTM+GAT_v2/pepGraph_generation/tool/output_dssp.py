#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os
from Bio import PDB
from Bio.PDB import PDBParser
from Bio.PDB import DSSP
from Bio.PDB.HSExposure import HSExposureCA


amino_acids = set([
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 
    'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 
    'SER', 'THR', 'TRP', 'TYR', 'VAL'
])
nucleotides = set(['DA', 'DC', 'DG', 'DT', 'DU', 'A', 'C', 'G', 'U'])

dssp_path = '/home/lwang/models/mambaforge/envs/liyao_env/bin/mkdssp'
def output_dssp(PDB_ID, PDB_file, save_dir, dssp_path = dssp_path):
    p = PDBParser()
    structure = p.get_structure(PDB_ID,PDB_file)
    model = structure[0]

    dssp = DSSP(model, PDB_file, dssp=dssp_path)
    chain_dssp = {}
    for key in dssp.keys():
        if key[0] not in chain_dssp.keys():
            chain_dssp[key[0]] = []
        chain_dssp[key[0]].append(dssp[key])
    
    for chain in chain_dssp.keys():
        with open(f'{save_dir}/{PDB_ID}_{chain}.dssp.txt',"w") as out:
            for aa_dssp in chain_dssp[chain]:
                for item in list(aa_dssp):
                    out.write(str(item) + "\t")
                out.write("\n")
    return chain_dssp


import biotite.structure as struc
import biotite.structure.io as strucio
import numpy as np
def biotite_SASA(pdb_file, chain_id, HDXparser = None):
    '''
    biotite.structure.sasa(array, probe_radius=1.4, atom_filter=None, ignore_ions=True, point_number=1000, point_distr='Fibonacci', vdw_radii='ProtOr')
    '''
    model = struc.io.load_structure(pdb_file)
    chain_mask = (model.chain_id == chain_id)
    backbone_mask = struc.filter_backbone(model)
    backbone_mask &= chain_mask

    atom_sasa = struc.sasa(model, atom_filter = backbone_mask, vdw_radii="Single")
    atom_sasa = np.nan_to_num(atom_sasa, nan = 0.0)

    return struc.apply_residue_wise(model, atom_sasa, np.sum), chain_mask

# from DeepRank
def get_bio_model(pdbfile):
    """Get the model

    Args:
        pdbfile (str): pdbfile

    Returns:
        [type]: Bio object
    """
    parser = PDBParser()
    structure = parser.get_structure('_tmp', pdbfile)
    return structure[0]

def get_hse(pdb_file, HDXparser):
    """Get the hydrogen surface exposure

    Args:
        model (bio model): model of the strucrture

    Returns:
        dict: hse data
    """
    model = get_bio_model(pdb_file)
    hse = HSExposureCA(model)
    data = {}
    hse_mtx = []
    index_range = []
    for k in list(hse.keys()):
        new_key = (k[0], k[1][1])
        index_range.append(k[1][1])

        x = hse[k]
        if x[2] is None:
            x = list(x)
            x[2] = 0.0
            x = tuple(x)

        data[new_key] = x
        hse_mtx.append(list(x))
    return data

def get_model_chains(structure):
    for model in structure:
        for chain in model:
            protein_residues = []
            nucleic_acid_residues = []
            ligands = []
            
            for residue in chain:
                # Use .get_resname() to get the residue name
                resname = residue.get_resname()
                
                if resname in amino_acids:
                    protein_residues.append(resname)
                elif resname in nucleotides:
                    nucleic_acid_residues.append(resname)
                else:
                    # Hetero residues not in standard amino acids or nucleotides are considered ligands
                    if residue.id[0] != ' ':
                        ligands.append(resname)
            
            if protein_residues:
                print(f"Chain {chain.id} is a protein chain with residues: {protein_residues[:5]}...")
            if nucleic_acid_residues:
                print(f"Chain {chain.id} is a nucleic acid chain with residues: {nucleic_acid_residues[:5]}...")
            if ligands:
                print(f"Chain {chain.id} contains ligands: {ligands}")

def parse_hhm(hhm_path):
    pass