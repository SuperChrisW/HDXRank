from MSA_embedding import generate_embedding
import os

root_dir = '/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/complex_model/'
dssp_list = ['LptDE_Apo.dssp.txt', 'LptDE_Thanatin.dssp.txt']
hhm_file = 'LptD.hhm'

os.chdir(root_dir)
generate_embedding(hhm_file, dssp_list[0], 'LptDE_Apo.embedding.txt')
generate_embedding(hhm_file, dssp_list[1], 'LptDE_Thanatin.embedding.txt')