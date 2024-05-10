#!/usr/bin/env python
__description__ = \
"""
preprocess.py
Find the match between atom numbers from two PDB files
"""
__author__ = "Kirill Speranskiy"
__date__ = "07-28-10"

import os
from pdb_data.common import *
import encode_hybrid36

# Path to instalation directory
homepath= ""

def setEnvVariable():
	"""
	Set environment variable where the newsim installed
	"""
	global homepath
	#for param in os.environ.keys():
	try: 
		homepath= os.environ["FRODANHOME"] +"/bin/"
		return
	except:
		print "\nEnvironment variable FRODANHOME is not set!"
		sys.exit()

def pdbSeq(pdb):
	"""
	Parse the ATOM entries in a pdb file.  Return dictionary of sequences keyed to chain.
	"""
	atoms=[]
        atoms_all = [l for l in pdb if l[0:6] == "ATOM  " or l[0:6] == "HETATM"]
	resnum_old = ""
	chain_id_old = ""
	for j in range(len(atoms_all) ): # Loop over PDB records
		atomline = atoms_all[j]
		resnum_new = atomline[23:26]
		if ( atomline[20] == " " ):
			chain_id_new= atomline[21]
		else:
			chain_id_new= atomline[20:22]
		if resnum_new != resnum_old or chain_id_new != chain_id_old:
			resnum_old = resnum_new
			chain_id_old = chain_id_new
			atoms.append(atomline)


        #chain_dict = dict([(l[21],[]) for l in atoms])
	chain_dict = dict()

	for atomline in atoms:
		if ( atomline[20] == " " ):
			chain_dict[ atomline[21] ]= []
		else:
			chain_dict[ atomline[20:22] ] = []

	for c in chain_dict.keys():
		#chain_dict[c] = [l[17:20] for l in atoms if l[21] == c]
		for atomline in atoms:
			if ( atomline[20] == " " and  atomline[21] == c):
				chain_dict[c].append( atomline[17:20] )
			elif ( atomline[20] != " " and  atomline[20:22] == c ):
				chain_dict[c].append( atomline[17:20] )
	return chain_dict

def pdbSeq2Fasta(pdb, pdb_id=""):
	"""
	Extract sequence from pdb file and write out in FASTA format.
	"""

	# Grab sequences
	chain_dict = pdbSeq(pdb)

	# Determine which chains are being written out
	chains_to_write = chain_dict.keys()

	# Convert sequences to 1-letter format and join strings
	for c in chains_to_write:
		for aa_index, aa in enumerate(chain_dict[c]):
			try:
				chain_dict[c][aa_index] = AA3_TO_AA1[aa]
			except KeyError:
				try: # If we cannot find amino acid try to find DNA
					chain_dict[c][aa_index] = DNA3_TO_DNA1[aa]
				except KeyError:
					chain_dict[c][aa_index] = "X"

	out = []
	out.append("> %s\n" % (pdb_id))
	for c in chains_to_write:
		# Write output in lines 80 characters long
		seq_length = len(chain_dict[c])
		num_lines = seq_length / 80
        
		for i in range(num_lines+1):
			out.append("".join([aa for aa in chain_dict[c][80*i:80*(i+1)]]))
			out.append("\n")
		#out.append("".join([aa for aa in chain_dict[c][80*(i+1):]]))
		#out.append("\n")

	return "".join(out)
 
           
def extract_sec(pdb_file1, pdb_file2):
	"""
	Function to extract sequence from PDB file.
	"""
        file_list=[]
        file_list.append(pdb_file1)
        file_list.append(pdb_file2)
        seqfile= "PDBaa.seq"
	fout = open(seqfile,'w')
	id_ind=1
	# Extract sequence data
	for pdb_file in file_list:
		pdb_id = os.path.split(pdb_file)[-1][:-4]
		pdb_id = "_%s_%s"%(pdb_id, id_ind)
                f = open(pdb_file,'r')
		pdb = f.readlines()
		f.close()
		seq = pdbSeq2Fasta(pdb, pdb_id)
		id_ind = id_ind+1
		print >>fout, seq
	fout.close()

   
# Parse the ATOM and HETATM records in PDB file
# Return dictionaries of atomic names and numbers keyed to serial number of residues
#
def Seq2Atom(pdb):
#	print pdb
	res_atomname_dict={} # Set correspondence of Residue serial number and atom names in this residue
	res_atomnumber_dict={} # Set correspondence of Residue serial number and atom serial number 
	atomname_array=[]
	atomnumber_array=[]
	ires_count=0
	atoms = [l for l in pdb if l[0:6] == "ATOM  " or l[0:6] == "HETATM"]
	atomline = atoms[0] # Read first line
	resnum_old = atomline[23:26] # Set first residue number
	chain_id_old = atomline[21]
	for j in range(len(atoms) ): # Loop over PDB records
		atomline = atoms[j]
		resnum_new = atomline[23:26]
		chain_id_new= atomline[21]
		if resnum_new != resnum_old or chain_id_new != chain_id_old:
			resnum_old = resnum_new
			chain_id_old = chain_id_new
			ires_count=ires_count+1
			atomname_array=[]
			atomnumber_array=[]
		atomname= atomline[12:16]
#		atomnum= int(atomline[6:11])-1
		atomnum=encode_hybrid36.hy36decode(5, atomline[6:11]) -1
		atomname_array.append(atomname)# Add atom name
		atomnumber_array.append(atomnum) # Add atom serial number 
		res_atomname_dict[ires_count] = atomname_array
		res_atomnumber_dict[ires_count] = atomnumber_array

	return res_atomname_dict, res_atomnumber_dict
	
# Find matching atoms in two structures using sequence alingment.
def MatchAtomNumber(file_list, seq_list):
	atomsMatch=[]
	res_atomname_dict={}
	res_atomnumber_dict={}
	# Extract atom data
	for pdb_file in file_list:
		pdb_id = os.path.split(pdb_file)[-1][:-4]
		f = open(pdb_file,'r')
		pdb = f.readlines()
		f.close()
		# List of atom names and serial numbers in each residue
		res_atomname_dict[pdb_file],res_atomnumber_dict[pdb_file]= Seq2Atom(pdb) 

	res_atomname_dict1= res_atomname_dict[file_list[0]] # Get dictionary of residue atom names from sequence1 
	res_atomnumber_dict1= res_atomnumber_dict[file_list[0]]# Get dictionary of residue atom numbers from sequence1 
	res_atomname_dict2= res_atomname_dict[file_list[1]]# Get dictionary of residue atom names from sequence2
	res_atomnumber_dict2= res_atomnumber_dict[file_list[1]]# Get dictionary of residue atom numbers from sequence2

	seq1= seq_list[0] # Array of one letter residues from sequence1
	seq2= seq_list[1] # Array of one letter residues from sequence2
	global iseq_len1
	global iseq_len2
	global	iinsertions
	global idelitions
	icount_first = 0
	icount_second = 0
	# Loop over aligned sequences
	for i in range ( len(seq1) ): # Length of seq1 and seq2 is the same at this point (they are aligned)
		atoms_res1=[] # Array of residue atoms from sequence1
		atoms_res2=[] # Array of residue atoms from sequence2
		res1= seq1[i]
		res2= seq2[i]
		control_flag=0
		if res1 != "-" and res2 != "-":
			atoms_res1.append( res_atomname_dict1[icount_first] )
			atoms_res1.append( res_atomnumber_dict1[icount_first] )
			atoms_res2.append( res_atomname_dict2[icount_second] )
			atoms_res2.append( res_atomnumber_dict2[icount_second] )
			icount_first = icount_first+1
			icount_second = icount_second+1
#			print "Case 1",  icount_first, icount_second, res1, res2
			control_flag=1 
		elif res1 != "-" and res2 == "-":
			icount_first = icount_first+1
#			print "Case 2", icount_first, icount_second, res1, res2			
		elif res1 == "-" and res2 != "-":
			icount_second = icount_second+1
#			print "Case 3",  icount_first, icount_second, res1, res2
		elif res1 == "-" and res2 == "-":
			print "Residues %s - %s are not defined in alignment file"%(res1, res2)
#			print "Case 4",  icount_first, icount_second, res1, res2			
		else:
			print "Error in sequence alignment"
		
		if control_flag:
			atoms1= atoms_res1[0]
			atoms2= atoms_res2[0]
			serial_atoms1= atoms_res1[1]
			serial_atoms2= atoms_res2[1]
			for iatom in range ( len(atoms_res1[0]) ):
				for jatom in range ( len(atoms_res2[0]) ):
					if atoms1[iatom] == atoms2[jatom]: # Compare atom names
						#print atoms1[iatom], serial_atoms1[iatom], atoms2[jatom], serial_atoms2[jatom]
						atomsMatch.append("%s %s"%(serial_atoms1[iatom], serial_atoms2[jatom]) )

	return atomsMatch

def split_pdb(pdb_file, ind):
	"""
	Split each PDB file in number of PDBs that contain only 1 chain. Return chain file names (_1_A.pdb, _1_B.pdb, _2_A.pdb, _2_B.pdb, etc. )
	"""
	
	chain_file_names=[]
        f = open(pdb_file,'r')
	pdb = f.readlines()
	f.close()
	atoms = [l for l in pdb if l[0:6] == "ATOM  " or l[0:6] == "HETATM" or l[0:3] == "TER"]
#	atoms_hash = [(l[0:6],i) for i, l in enumerate(atoms)]
        atoms_hash = []
        for i, l in enumerate(atoms):
            if(l[0:6] == "ATOM  " or l[0:6] == "HETATM"):
                atoms_hash.append((l[0:6],i))
            if (l[0:3] == "TER" ):
                atoms_hash.append((l[0:3],i))

        atoms_hash_1 = [x[1] for x in atoms_hash if x[0] == "ATOM  " or x[0] == "HETATM"]         
	atoms_hash_2 = [x[1] for x in atoms_hash if x[0] == "TER"]
	last_record=atoms[ len(atoms_hash)-1]
	if (last_record[0:3] != "TER"):
		atoms_hash_2.append(len(atoms_hash))
	atoms_hash_2.append(-1)
	
	all_chains = []
        all_chains.append(pdb[atoms_hash_1[0]:atoms_hash_2[0]])
	for i in range(1,len(atoms_hash_2)-1):
		all_chains.append(pdb[atoms_hash_2[i-1]:atoms_hash_2[i]])
	
        for index, chain in enumerate(all_chains):
		chain_dict = pdbSeq(chain)
		chain_key = chain_dict.keys()
		filename= "_%i_%s.pdb" % (ind, chain_key[0])
		g = open(filename,"w")
		chain_file_names.append( filename )
		g.writelines(chain)
		g.close()
	
	return chain_file_names 

def parseAlignment():
	global homepath
	setEnvVariable()	
	command= homepath + "alignment_amino PDBaa.seq PDBaligned.seq"
	os.system(command)	
	f = open("PDBaligned.seq",'r')
	lines= f.readlines()
	f.close()
	icount=0
	seq1=""
	seq2=""
	for line in lines:
		if (len(line) >1 ):
			if (icount == 1 ):
				seq1 +=line.strip()
			if (icount == 3 ):
				seq2 +=line.strip()
			icount +=1
		else:
			icount =0
	seq_list=[]
	seq_list.append(seq1)
	seq_list.append(seq2)

	return seq_list, lines
	
def mapAtoms(file_list, chain_file_name):
	"""
	Function to execute if called from preprocess.py. 
	"""
	chainfile= chain_file_name
                                                                               
	pdb_chains = {}
	ind_pdb =1
	for pdb_file in file_list: #Split each PDB file in number of PDBs that contain only 1 chain (with names _1_A, _1_B,etc. and _2_A, _2_B, etc.)
		pdb_chains[ind_pdb] = split_pdb(pdb_file, ind_pdb)
                ind_pdb = ind_pdb+1

        f = open(chainfile,'r')
	chains = f.readlines()
	f.close()

	chain_init=[]

	chain_init= chains[0].split()
	chain_targ=[]
	chain_targ= chains[1].split()
        atomMatch=[]
	alinment=[] # List that will contain alinment of all chains
	str=''
        for i in range (len(chain_init)):
                if (chain_init[i] != "-" and chain_targ[i] != "-"):
		  str= "Chains: %s %s \n" % (chain_init[i], chain_targ[i])
                  file_name1="_1_%s.pdb" % (chain_init[i])
                  file_name2="_2_%s.pdb" % (chain_targ[i])
	          extract_sec( file_name1, file_name2 )
		  seq_list,  chain_alinment= parseAlignment()
                  file_list_chains=[]
	          file_list_chains.append(file_name1)
                  file_list_chains.append(file_name2)
	          atomMatch_chain= MatchAtomNumber(file_list_chains, seq_list)
                  atomMatch.append(atomMatch_chain)
		  alinment.append(str)
		  for i in range ( len(chain_alinment)):
			  alinment.append(chain_alinment[i])
		  
	fout = open("PDBaligned.seq",'w')
	for line in alinment:
		print>>fout, line,
	fout.close()	
	out=[]
	    
        for i, aa in enumerate( pdb_chains ):	
		for ii in range ( len(pdb_chains[aa])):
			os.system("rm -f "+ pdb_chains[aa][ii])

	atom_matchfile ="atommatch.txt"	
	for i in range( len(atomMatch) ):
                atomMatch_chain= atomMatch[i]
                for j in range ( len(atomMatch_chain) ):
		  out.append(atomMatch_chain[j])
		  out.append("\n")
	
        fout = open(atom_matchfile,'w')
	print >>fout, "".join(out)	
	fout.close()
	
	print "Finished mapping atoms. Result file:",atom_matchfile
	return atom_matchfile
