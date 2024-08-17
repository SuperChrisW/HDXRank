import os
import subprocess
import shutil
import xml.etree.ElementTree as ET

# automatic run first analysing proflex dataset, output rigid cluster file

def read_RC(key, direct):
    clst = {}
    path = os.path.join(direct,key)

    with open(path, 'r') as f:
        data = f.read().strip().split('\n') 
        for line in data:
            numbers = line.strip().split()
            atom_idx, aname, cluster = numbers[0].strip(), numbers[1].strip(), numbers[2].strip()
            if "X" not in aname:
                if cluster not in clst.keys():
                    clst[cluster] = []
                clst[cluster].append(atom_idx)
    return clst

def RUN_FIRST(proflex_dir, key, para):
    filename = key
#    filepath = os.path.join(dict, filename)
#    print(filepath)

    # Set PROFLEX_HOME
    os.chdir(proflex_dir)
    os.environ["PROFLEX_HOME"] = f"{proflex_dir}/prog"
    os.environ["PATH"] += os.pathsep + f"{proflex_dir}/bin"

    if os.path.isfile(filename):
        command = f'proflex {para} {filename}'
        # conduct the terminal command for each protein in temp_list
        subprocess.run(command, shell=True)
    else:
        print(f'{filename} does not exist')
    return None

def RUN_FRODAN(key, para):
    filename= key.split('.') # assume key contains suffix
    PDB = filename[0]
    os.chdir('/home/wang/Documents/Frodan/frodaN')
    os.environ["FRODANHOME"] = "/home/wang/Documents/Frodan/frodaN"

    if os.path.isfile(filename):
        tree = ET.parse('options_fixedcons.xml')
        root = tree.getroot()

        # Find the 'modelfiles' element
        modelfiles = root.find('./modelfiles')
        if modelfiles is not None:
            # Set new attribute values
            modelfiles.set('pdb', key)
            modelfiles.set('cov', f'cov_{PDB}.txt')
            modelfiles.set('rigidunits', f'rc_{PDB}.txt')
            modelfiles.set('atomtypes', f'atomtypes_{PDB}.txt')

            option_file = f'options_{PDB}.xml'
            tree.write(option_file)

        else:
            print("The <modelfiles> tag was not found!")
            return False

        command = f'python2 SFRODANHOME/preprocess.py -i {key} -o options_{PDB}.xml'
        # conduct the terminal command for each protein in temp_list
        subprocess.run(command, shell=True)

def list_files_in_directory(directory):
    return [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]

def copyfiles(key, proflex_dir, dst_dir): 
    proflex_key = f'{key}_proflexdataset'
    cluster_key = f'decomp_list'
    src_cluster = os.path.join(proflex_dir, cluster_key)
    src_proflex = os.path.join(proflex_dir, proflex_key)
    
    path1 = f'/decomp_{key}'
    path2 = f'/{proflex_key}'
    cluster_filepath = dst_dir + path1
    proflex_filepath = dst_dir + path2
    
    if os.path.isfile(src_proflex) and os.path.isfile(src_cluster):
        shutil.move(src_cluster, cluster_filepath)        
        shutil.move(src_proflex, proflex_filepath)
        return True
    else:
        return False
  

def pdb_clean(src_path, dst_path):
    if os.path.exists(src_path):
        with open(src_path, 'r') as src_file, open(dst_path, 'w') as dst_file:
            dst_file.write('HEADER    \n')
            for line in src_file:
                if line.startswith('ATOM') or line.startswith('HETATM') or line.startswith('TER'):
                    dst_file.write(line)
            dst_file.write('END')
        return True
    else:
        print("Source file does not exist.")
        return False
    
def hbplus(hbplus_dir, filename, save_dir):
    os.environ['PATH'] += os.pathsep + hbplus_dir
    os.chdir(hbplus_dir)

    filename = filename.split('.')[0]
    src_file = os.path.join(hbplus_dir, f'{filename}.pdb')
    command = f'hbplus -O {src_file}'
    os.system(command)
    os.remove(src_file)

    src_file = os.path.join(hbplus_dir, f'{filename}.h')
    dst_file = os.path.join(save_dir, f'{filename}_Hplus.pdb')
    shutil.move(src_file, dst_file)

    src_file = os.path.join(hbplus_dir, f'{filename}.hb2')
    dst_file = os.path.join(save_dir, f'{filename}.hb2')
    shutil.move(src_file, dst_file)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           