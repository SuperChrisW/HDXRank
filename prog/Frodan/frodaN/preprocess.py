#!/usr/bin/env python
__description__ = \
"""
preprocess.py
Prepare input files for GEMS
"""
__author__ = "Kirill Speranskiy"
__date__ = "05-08-2010"

import os
import sys
import xml.dom.minidom
from xml.dom.minidom import Node
import  python_scripts.map_atoms
import  python_scripts.validate_structure
import python_scripts.encode_hybrid36


# Path to instalation directory
homepath= ""

# 
def checkRequiredFiles(files): 
	"""
	Files tag in options.xml should have all atributes. These atributes specify all files used by gems.
	"""
	filenames=["pdb", "cov", "rigidunits", "atomtypes"]
	for filename in filenames:
		if filename not in files:
			print filename + " file is not found. Please add "+filename + " to the list of files"
			print "Current list of files:  "
			for s, ss in enumerate(files):
				print "%s = \"%s\" "%(ss, files[ss])
			sys.exit()
	print "All file names are found. Starting external binaries ..."

		
def setEnvVariable(): 
	"""
	Set environment variable where the gems installed
	"""
	global homepath
	#for param in os.environ.keys():
	try: 
		homepath= os.environ["FRODANHOME"]+"/bin/"
		return
	except:
		print "\nEnvironment variable FRODANHOME is not set!"
		sys.exit()


def parseParamenters(xml_file, tagname): 
	"""
	Parse parameter xml file of the run
	"""
	tagname_map={}
	doc = xml.dom.minidom.parse(xml_file)
	#optionsNode= doc.firstChild
	#for child in optionsNode.childNodes:
	#	if child.nodeType == child.ELEMENT_NODE :
	#		if child.tagName == tagname:
	for node in doc.getElementsByTagName(tagname):
					attrs= node.attributes
					for i in range (attrs.length):
						attr= attrs.item(i)
						tagname_map[str(attr.name)] = str(attr.value)
			
	
	checkRequiredFiles(tagname_map)
	return tagname_map


def splitName(pdbname):
	"""
	Splits PDB name
	"""
	name= pdbname.split(".")
	basename= name[0]
	ext= name[1]
	
	return basename, ext

def runFirst(pdbname, files_dict): # Run FIRST and move resulting files to specified names
	global homepath
	basename, ext = splitName(pdbname)
	command = homepath+"FIRST "+ pdbname + " -non  -numbers -nohbonds -nohphobes -noaromatics -covout -L %s/../FIRST/"%(homepath)
	os.system(command)
	basename, ext= splitName(pdbname)
	command = "mv %s_RCD.pdb %s"%(basename, files_dict["pdb"])
	os.system(command)
	command = "mv cov.out %s"%(files_dict["cov"])
	os.system(command)
	command = "mv %s_data.txt %s"%(basename, files_dict["rigidunits"])
	os.system(command)
		
def runAtomTypes(file_dict):
	global homepath
	command = homepath+"assignAtomTypes --pdb %s --cov %s -o %s --FIRSTBondFileInterpretation index1"%(file_dict["pdb"], file_dict["cov"],file_dict["atomtypes"] )
	os.system(command)

def runMatchState(init_files, targ_files):
	global homepath	
	#Run matchState1ToState2
	command = homepath+"matchState1ToState2 --pdb1 %s --FIRSTcov1 %s --FIRSTFileFormat index1 --pdb2 %s -m %s --out1 %s"%(init_files["pdb"], init_files["cov"], targ_files["pdb"], targ_files["map"], init_files["pdb"])
	os.system(command)	

def check_files(file_name_init, hetatm): # Run preliminary check of PDB file format
	basename, ext = splitName(file_name_init)
       	file_name = basename+"_h36."+ext 
	f = open(file_name_init,'r')
	pdb = f.readlines()
	pdb_validated = python_scripts.validate_structure.validate(pdb, hetatm)
	pdb_encoded= python_scripts.encode_hybrid36.encode36(pdb_validated)
	f.close()
	f = open(file_name,'w')
	for pdb_line in pdb_encoded:
		print >> f, pdb_line,
	f.close()
	
	return file_name

def parseCommandLine(argv):
	comline_arg={}
	for s in argv:
		if s ==  "-i":
			ind = argv.index(s)
			comline_arg[s] = argv[ind+1]
			if (os.path.isfile(argv[ind+1]) != True):
				print "File "+ argv[ind+1] +" not found."
				usage()
	
		if s ==  "-t":
			ind = argv.index(s)
			comline_arg[s] = argv[ind+1]
			if (os.path.isfile(argv[ind+1]) != True):
				print "File "+ argv[ind+1] +" not found."
				usage()
	
		if s ==  "-o":
			ind = argv.index(s)
			comline_arg[s] = argv[ind+1]
			if (os.path.isfile(argv[ind+1]) != True):
				print "File "+ argv[ind+1] +" not found."
				usage()

		if s ==  "-nohet":
			ind = argv.index(s)
			comline_arg[s] = argv[ind]
			
	if "-i" not in comline_arg:	
		usage()
	if "-o" not in comline_arg:
		usage()
	if "-help" in comline_arg or ( len(argv) == 1):
		usage()

	return comline_arg

def usage():
	print "Usage: %s  -i initial.pdb [-t final.pdb] -o options.xml [-nohet]"%(sys.argv[0])
	print "-i set initial PDB structure"
	print "-t set target PDB structure. This option is required for targeting run"
	print "-o set run options XML file"
	print "-nohet do not include heterogens"
	print "-help print this help"
	sys.exit()
		
def cleanOutput(file_list, comline_arg):
	for file in file_list:
		basename, ext = splitName(file)
		command="rm %s* "%(basename)
		os.system(command)
	if "-t" in comline_arg:
		command="rm PDBaa.seq"
		os.system(command)
	
def main():
	"""
	Function to execute if called from command line. 
	Command: preprocess.py -i initial.pdb -t final.pdb -o options.xml -nohet
	"""
	setEnvVariable()
	comline_arg= parseCommandLine(sys.argv) # Parse command line arguments

	hetatm= True # By default HETATM atoms are included in calculations
	if "-nohet" in comline_arg:
		hetatm= False
			
	file_list=[]
	if "-i" in comline_arg:
		file_list.append("")
		init_files = parseParamenters(comline_arg["-o"], "modelfiles")
	        filename = check_files(comline_arg["-i"], hetatm)
		file_list[0]= filename 
		if (comline_arg["-i"] == init_files["pdb"]):
			print "Initial and processed file names are the same. Please change processed file name."
			sys.exit()
		runFirst(filename, init_files)
		runAtomTypes(init_files)
		
	if "-t" in comline_arg:
		file_list.append("")
		targ_files = parseParamenters(comline_arg["-o"], "targetfiles")
	        filename = check_files(comline_arg["-t"], hetatm)
		file_list[1]=filename 
		if (comline_arg["-t"] == init_files["pdb"]):
			print "Target and processed file names are the same. Please change processed file name."
			sys.exit()
		runFirst(filename, targ_files)
		runAtomTypes(targ_files)
		if "map" in targ_files:
			atom_matchfile = python_scripts.map_atoms.mapAtoms(file_list)
			command = "mv %s %s"%(atom_matchfile, targ_files["map"])
			os.system(command)
			runMatchState(init_files, targ_files)
		else:
			print "Map file is not specified in targetfiles field."	

	cleanOutput(file_list, comline_arg)

	
if __name__ == "__main__":
    main()
