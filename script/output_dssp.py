#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 15:38:40 2021

@author: jiali
Calculate surface accessibility from PDB
Usage: python output_dssp.py <ID> <PDB_file.pdb>
"""
import sys
import os
from Bio import PDB
from Bio.PDB import PDBParser
from Bio.PDB import DSSP
from predataprocess import get_pdb_seq

dssp_path = '/home/lwang/miniconda3/envs/liyao_env/bin/mkdssp'

def output_dssp(PDB_ID, PDB_file, chain_id, save_dir):
    p = PDBParser()
    structure = p.get_structure(PDB_ID,PDB_file)
    model = structure[0]
    # calculate the number of AA in the protein
    aa_list = []
    dssp = DSSP(model, PDB_file, dssp=dssp_path)

    '''
    res_start, res_end = 0, 0
    pos_record = {}
    for chain in list(structure.get_chains()):
        sequence, res_list = get_pdb_seq(structure, chain.id, fill_gap=False)
        res_end = res_start + len(sequence) # FIXME: this is not correct
        pos_record[chain.id] = (res_start, res_end)
        res_start = res_end
        for id in chain.get_residues():
            if PDB.is_aa(id):
                aa_list.append(id)
    print(pos_record)

    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    start_idx = pos_record[chain_id][0]
    end_idx = pos_record[chain_id][1]
    '''

    chain_dssp = []
    for key in dssp.keys():
        if key[0] == chain_id:
            chain_dssp.append(dssp[key])
    
    with open(f'{save_dir}/{PDB_ID}.dssp.txt',"w") as out:
        for aa_dssp in chain_dssp:
            for item in list(aa_dssp):
                out.write(str(item) + "\t")
            out.write("\n")

    '''
    with open(f'{save_dir}/{PDB_ID}.dssp.txt',"w") as out:
        start_idx = pos_record[chain_id][0]
        end_idx = pos_record[chain_id][1]
        
        for a_key in list(dssp.keys())[start_idx:end_idx]:
            for item in list(dssp[a_key]):
                out.write(str(item) + "\t")
            out.write("\n")
    '''