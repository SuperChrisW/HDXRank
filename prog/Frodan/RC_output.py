import subprocess
import shutil
import xml.etree.ElementTree as ET
import os

def RUN_FRODAN(key):
    filename = key.split('.')  # assume key contains suffix
    PDB = filename[0]
    os.chdir('/home/wang/Documents/Frodan/frodaN')
    os.environ["FRODANHOME"] = "/home/wang/Documents/Frodan/frodaN"

    if os.path.isfile(key):
        tree = ET.parse('options_fixedcons.xml')
        root = tree.getroot()

        # Find the 'modelfiles' element
        modelfiles = root.find('./modelfiles')
        if modelfiles is not None:
            # Set new attribute values
            modelfiles.set('pdb', '{}_process.pdb'.format(PDB))
            modelfiles.set('cov', 'cov_{}.txt'.format(PDB))
            modelfiles.set('rigidunits', 'rc_{}.txt'.format(PDB))
            modelfiles.set('atomtypes', 'atomtypes_{}.txt'.format(PDB))

            option_file = 'options_{}.xml'.format(PDB)
            tree.write(option_file)

        else:
            print "The <modelfiles> tag was not found!"
            return False
         
    command1 = 'python2 $FRODANHOME/preprocess.py -i {} -o options_{}.xml'.format(key, PDB)
    command2 = '$FRODANHOME/bin/main options_{}.xml'.format(PDB)
    subprocess.call(command1, shell=True)  # changed from subprocess.run to subprocess.call
    subprocess.call(command2, shell=True)
    return True

def list_files_in_directory(directory):
    return [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]

def copyfiles(key, proflex_dir, dst_dir): 
    if not os.path.exists(dst_dir):
    	os.mkdir(dst_dir)
    else:
    	return False
    filename = key.split('.')  # assume key contains suffix
    PDB = filename[0]
    
    file_list = ['rc_{}.txt'.format(PDB), 'options_{}.xml'.format(PDB), 'cov_{}.txt'.format(PDB), 'atomtypes_{}.txt'.format(PDB)]
    for file in file_list:
	src_cluster = os.path.join(proflex_dir, file)
	cluster_filepath = os.path.join(dst_dir, file)
        shutil.move(src_cluster, cluster_filepath)        

    command = 'mv {}/{}*.pdb {}/'.format(proflex_dir, PDB, dst_dir)
    subprocess.call(command, shell = True)
    return True
