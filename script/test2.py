## output residue ANSSURR rigidity ##
from prot_rigidity import ASSURR_flexibility

root_dir = '/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/HDX_MS_dataset/AF_model'

ASSURR_flexibility(f'{root_dir}/AF_proflex', f'{root_dir}/AF_dssp', f'{root_dir}/AF_rigidity')