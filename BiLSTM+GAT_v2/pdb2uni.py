from urllib.request import urlopen
from Bio import SwissProt
import json
import pandas as pd

def get_uniprot_sequence(uniprot_id):
    try:
        record = SwissProt.read(SwissProt.urlopen(f"https://www.uniprot.org/uniprot/{uniprot_id}.txt"))
        return record.sequence
    except Exception as e:
        print(f"Error fetching data for UniProt ID {uniprot_id}: {str(e)}")
        return None

# read xlsx file
df = pd.read_excel('/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/HDX_MS_dataset/dataset_summary.xlsx', sheet_name='Sheet3')
pdb_col = 'PDB_ID'
chain_col = 'Chain'
pdb_chain_df = df[[pdb_col, chain_col]]

if 'uniprot' in pdb_chain_df.columns:
    uniprot_ids = pdb_chain_df['uniprot']
else:
    uniprot_ids = []

result = {'pdb':[], 'chain':[], 'uniprot':[], 'sequence':[]}

for pdb, chain, uniprot_id in zip(pdb_chain_df[pdb_col], pdb_chain_df[chain_col], uniprot_ids):
    print('mapping...', pdb, chain)
    try:
        content = urlopen('https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/' + pdb).read()
    except:
        print(pdb, chain, "PDB Not Found (HTTP Error 404). Skipped.")
        continue
    
    content = json.loads(content.decode('utf-8'))

    # check if chain exists
    chain_exist_boo_list = []

    # find uniprot id
    for uniprot in content[pdb.lower()]['UniProt'].keys():
        for mapping in content[pdb.lower()]['UniProt'][uniprot]['mappings']:
            if mapping['chain_id'] == chain.upper():
                uniprot_id = uniprot if uniprot else uniprot_id
                uniprot_sequence = get_uniprot_sequence(uniprot)

                if uniprot_sequence:
                    result['pdb'].append(pdb)
                    result['chain'].append(chain)
                    result['uniprot'].append(uniprot_id)
                    result['sequence'].append(uniprot_sequence)
                chain_exist_boo_list.append(True)
            else:
                chain_exist_boo_list.append(False)

    if not any(chain_exist_boo_list):
        print(pdb, chain, "PDB Found but Chain Not Found. Skipped.")

result_df = pd.DataFrame(result)

csv_df = df.set_index([pdb_col, chain_col])
result_df.columns = [pdb_col, chain_col, 'Uniprot_ID', 'Sequence']
result_df = result_df.set_index([pdb_col, chain_col])
csv_df = csv_df.join(result_df)
# Save the merged DataFrame to a new CSV file
csv_df.to_csv('HDX_uniprot.csv', na_rep='NaN')