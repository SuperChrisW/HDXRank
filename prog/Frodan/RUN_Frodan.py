import os
import subprocess
import shutil
from RC_output import *

wild_path = '/home/wang/Documents/Archive'
#mut_path = '/home/wang/Documents/mut_pdb_addH'
frodan_path = '/home/wang/Documents/Frodan/frodaN'
rigid_path = '/home/wang/Documents/HDX_Dynamic'

wild_list = list_files_in_directory(wild_path)
#mut_list = list_files_in_directory(mut_path)
pdb_lines = []
count_wild, count_mut = 0, 0

for file in wild_list:
    key, suffix1, suffix2 = file.split('.')
    wildsave_path = os.path.join(rigid_path, key)  # changed from f-string to os.path.join
    if os.path.exists(wildsave_path): continue
    
    rename = key + '.pdb'
    src_path = os.path.join(wild_path, file)
    dst_path = os.path.join(frodan_path, rename)

    shutil.copy2(src_path, dst_path)
    RUN_FRODAN(rename)
    if copyfiles(rename, frodan_path, wildsave_path):
        count_wild += 1


'''
for file in mut_list:
    key, suffix = file.split('.')
    src_path = os.path.join(mut_path, file)
    dst_path = os.path.join(frodan_path, file)

    shutil.copy2(src_path, dst_path)
    RUN_FRODAN(file)
    mutsave_path = os.path.join(rigid_path, 'mut')  # changed from f-string to os.path.join
    if copyfiles(key, frodan_path, mutsave_path):
        count_mut += 1

    pdb_lines = []  # changed from clear() method
'''

print count_wild
#print count_mut
