
"""
Created on Wed Jan. 24 2024
@author: Liyao
Assemble all features into one embedding file
"""
import torch
import os
from Bio import PDB
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Polypeptide import one_to_three as protein_letters_1to3
import pandas as pd
import numpy as np

from dataclasses import dataclass
from pepGraph_generation.pepGraph_utlis import generate_embedding, AA_array, residue_charge, get_seq_polarity, parse_hhm
from output_dssp import biotite_SASA, get_hse, get_bio_model, get_model_chains


@dataclass
class RawInputData:
    msa: torch.Tensor
    ins: torch.Tensor
    res_HDMD: torch.Tensor
    res_polarity: torch.Tensor
    res_charge: torch.Tensor
    dssp: torch.Tensor
    rigidity: torch.Tensor
    hse: torch.Tensor
    log_t: torch.Tensor
    seq_len: torch.Tensor

    def construct_embedding(self):
        # construct embedding
        embedding = torch.cat((self.msa, self.ins, self.res_HDMD, self.res_polarity, self.res_charge, self.dssp, self.rigidity, self.hse, self.log_t), dim=1)
        return embedding

class HDXparser:
    def __init__(self, HDX_datatable, root_dir):
        self.HDX_datatable = HDX_datatable
        self.root_dir = root_dir

        self.dssp_dir = os.path.join(root_dir, 'dssp_files')
        #self.DSSP = "/home/lwang/miniconda3/envs/liyao_env/bin/mkdssp"
        self.hhm_dir = os.path.join(root_dir, 'hhm_files')
        self.rigidity_dir = os.path.join(root_dir, 'rigidity_files')
        self.proflex_dir = os.path.join(root_dir, 'proflex_files')
        self.pdb_dir = os.path.join(root_dir, 'structure')
        self.save_dir = os.path.join(root_dir, 'embedding_files')

    def construct_feat(self):
        temp_HDX_datatable = self.HDX_datatable.drop_duplicates(subset=['apo_identifier'])
        for index, row in temp_HDX_datatable.iterrows():
            pdbfile = os.path.join(self.pdb_dir, f'{row["apo_identifier"].str.strip()}')
            Bio_model = get_bio_model(pdbfile)
            chain_dict = get_model_chains(Bio_model) # get all chains in the model including Nucleic acid and small molecules

            # SASA
            biotite_model, atom_sasa, res_sasa = biotite_SASA(pdbfile, self)
            chain_mask = biotite_model.chain_id 
            # HSE
            res_hse = get_hse(pdbfile, self)
            # Rigidity
            pass
            
            # Res HDMD/polarity/charge
            ppb = PPBuilder()
            pdb_seq = ''
            for pp in ppb.build_peptides(Bio_model, aa_only=False):
                pdb_seq += str(pp.get_sequence())
            res_charge = np.array([residue_charge[protein_letters_1to3(res)] for res in pdb_seq]).reshape(-1, 1)
            res_polarity = get_seq_polarity(pdb_seq)
            HDMD = np.array([AA_array[res] for res in pdb_seq]).reshape(-1, 5)

            # HHblits MSA
            hhm_mtx, hhm_seq = parse_hhm()

        return load_protein(
            sasa,
            res_hse,
            HDMD,
            res_polarity,
            res_charge,
            hhm_file
        )


def main():
    HDX_datatable = pd.read_csv('HDX_datatable.csv')
    root_dir = '/home/lwang/AI-HDX-main/ProteinComplex_HDX_prediction'
    HDXparser = HDXparser(HDX_datatable, root_dir)
    InputData = HDXparser.construct_feat()

    embedding = InputData.construct_embedding()
    torch.save(embedding, HDXparser.save_dir)
    print(f'Embedding file is saved at {HDXparser.save_dir}')

if __name__ == "__main__":
    main()