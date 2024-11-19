from pepGraph_dataset import *
from pepGraph_utlis import load_embedding
import torch_geometric as tg
from torch.utils.data import Dataset
import torch

class data_process(Dataset):
    def __init__(self, **kwargs):
        self.max_len = 30
        pass

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):
        pass

    @staticmethod
    def read_HDX_table(HDX_df, protein, state, chain,  vector_file = None, pdb_file = None,
                    mismatch_report = True, correction_value = 0, filter = []):
        '''
        read the HDX table file and return the residue information as node_list
        '''
        pep._registry = [] # reset the peptide registry
        npep = 0
        res_chainsplit = {}
        for residue in res._registry:
            if residue.chain_id not in res_chainsplit.keys():
                res_chainsplit[residue.chain_id] = []
            res_chainsplit[residue.chain_id].append(residue)

        ## seperately process states
        for index, (protein, state, chain) in enumerate(zip(proteins, states, chains)):
            temp_HDX_df = HDX_df[(HDX_df['state']==state) & (HDX_df['protein']==protein)]
            #temp_HDX_df = temp_HDX_df[temp_HDX_df['sequence'].str.len() < max_len]
            temp_HDX_df = temp_HDX_df.sort_values(by=['start', 'end'], ascending=[True, True])
            #temp_HDX_df = temp_HDX_df.drop_duplicates(subset=['sequence'])

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