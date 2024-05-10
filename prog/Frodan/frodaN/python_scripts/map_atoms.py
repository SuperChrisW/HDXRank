#!/usr/bin/env python
__description__ = \
"""
preprocess.py
Find the match between atom numbers from two PDB files
"""
__author__ = "Kirill Speranskiy"
__date__ = "04-28-10"


import os
from pdb_data.common import *

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
	chains1=[]
	chain_arr_keys=[] # Arrays contains the correct sequence of keys for chain_dict
	for j in range(len(atoms_all) ): # Loop over PDB records
		atomline = atoms_all[j]
		resnum_new = atomline[22:26]
	       	if ( atomline[20] == " " ):
			chain_id_new= atomline[21]
		else:
			chain_id_new= atomline[20:22]
			
		if chain_id_new != chain_id_old:
			chain_arr_keys.append(chain_id_new)
		if resnum_new != resnum_old or chain_id_new != chain_id_old:
			resnum_old = resnum_new
			chain_id_old = chain_id_new
			atoms.append(atomline)
			chains1.append(chain_id_new)
			#chain_set.add(chain_id_new)
		
	chain_dict = dict()

        for c in chain_arr_keys:
		chain_dict[c]= []

        for c in chain_arr_keys:
             for atomline in atoms:
			if ( atomline[20] == " " and  atomline[21] == c):
				chain_dict[c].append( atomline[17:20] )
			elif ( atomline[20] != " " and  atomline[20:22] == c ):
				chain_dict[c].append( atomline[17:20] )

	dict_count=dict()
	for c in chain_arr_keys:
		dict_count[c]=0

	return chain_dict, chain_arr_keys


def pdbNumbers(pdb):
	"""
	Parse the ATOM entries in a pdb file.  Return dictionary of residue numbers keyed to chain.
	"""
	atoms=[]
        atoms_all = [l for l in pdb if l[0:6] == "ATOM  " or l[0:6] == "HETATM"]
	resnum_old = ""
	chain_id_old = ""
	for j in range(len(atoms_all) ): # Loop over PDB records
		atomline = atoms_all[j]
		resnum_new = atomline[22:26]

		if ( atomline[20] == " " ):
			chain_id_new= atomline[21]
		else:
			chain_id_new= atomline[20:22]

		if resnum_new != resnum_old or chain_id_new != chain_id_old:
			resnum_old = resnum_new
			chain_id_old = chain_id_new
			atoms.append(atomline)

	numbers_dict = dict()
	for atomline in atoms:
		if ( atomline[20] == " " ):
			numbers_dict[ atomline[21] ]= []
		else:
			numbers_dict[ atomline[20:22] ] = []

        for c in numbers_dict.keys():
	    for atomline in atoms:
			if ( atomline[20] == " " and  atomline[21] == c):
				numbers_dict[c].append( atomline[22:26] )
			elif ( atomline[20] != " " and  atomline[20:22] == c ):
				numbers_dict[c].append( atomline[22:26] )

	return numbers_dict 

def pdbSeq2Fasta(pdb, pdb_id=""):
	"""
	Extract sequence from pdb file and write out in FASTA format.
	"""

	# Grab sequences
	chain_dict, chain_arr_keys = pdbSeq(pdb)
	# Determine which chains are being written out
	chains_to_write = chain_arr_keys

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

	return "".join(out)
 
           
def extract_seq(pdb_file1, pdb_file2):
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
		pdb_id = "%s_%s"%(pdb_id, id_ind)
                f = open(pdb_file,'r')
		pdb = f.readlines()
		f.close()
		seq = pdbSeq2Fasta(pdb, pdb_id)
		id_ind = id_ind+1
		print >>fout, seq
	fout.close()
	
        print "Finished extracting sequence"
   
# Try to align sequences using residue numbers from PDB file.
# If residue numbers match in PDB structures the residue IDs will also match, otherwise we will put "dash" in the sequences
# This function is available only if user uses ClustalW
def do_rough_alingment(pdb_file1, pdb_file2):

        file_list=[]
	keys_arr=[]
        file_list.append(pdb_file1)
        file_list.append(pdb_file2)
        seqfile= "PDBaa.gde"

        seq_array=[]
	numbers_array=[]
	fout = open(seqfile,'w')
	id_ind=1
	# Extract sequence data
	for pdb_file in file_list:
                f = open(pdb_file,'r')
		pdb = f.readlines()
		f.close()
		chain_dict,chain_arr_keys = pdbSeq(pdb)
		numbers_dict = pdbNumbers(pdb)
		#print numbers_dict
		id_ind = id_ind+1
		seq_array.append(chain_dict)
		numbers_array.append(numbers_dict)
		keys_arr.append(chain_arr_keys)

	for i in range ( len(seq_array) ):
		chain_dict = seq_array[i]
		# Determine which chains are being written out
		chains_to_write = chain_dict.keys()
		# Convert sequences to 1-letter format
		for c in chains_to_write:
			for aa_index, aa in enumerate(chain_dict[c]):
				try:
					chain_dict[c][aa_index] = AA3_TO_AA1[aa]
				except KeyError:
					try:
						chain_dict[c][aa_index] = DNA3_TO_DNA1[aa]
					except KeyError:
						chain_dict[c][aa_index]  = "X"


	chain_dict1 = seq_array[0]
	numbers_dict1 = numbers_array[0]
	numbers_dict1_new={}
	chains_array1_newkeys=[]
	chains_array2_newkeys=[]	
	
	chain_dict2 = seq_array[1]	
	numbers_dict2 = numbers_array[1]
	numbers_dict2_new={}
	chains_dict1_new={}
	chains_dict2_new={}
	# Determine which chains are being written out in first structure
	#chains_to_write1 = chain_dict1.keys()
	chains_to_write1 =keys_arr[0]
	# Convert array to dictionary
	count_chains=0
	for c in chains_to_write1:
		for nn_index, nn in enumerate( numbers_dict1[c] ):
			numbers_dict1_new["%s,%s"%(c,nn)] = numbers_dict1[c][nn_index]
			chains_dict1_new["%s,%s"%(c,nn)]=chain_dict1[c][nn_index]
			chains_array1_newkeys.append("%s,%s"%(c,nn) )
		count_chains +=1

	chains_to_write2 =keys_arr[1]

	# Convert array to dictionary 
	for c in chains_to_write2:		
		for nn_index, nn in enumerate( numbers_dict2[c] ):
			numbers_dict2_new["%s,%s"%(c,nn)] = numbers_dict2[c][nn_index]
			chains_dict2_new["%s,%s"%(c,nn)]=chain_dict2[c][nn_index]
			chains_array2_newkeys.append("%s,%s"%(c,nn) )
			
	# Add missing elements in each dictionary
	dict1_keys= numbers_dict1_new.keys()
	i=0
	for kk in chains_array1_newkeys:
		try:
			numbers_dict2_new[kk]
		except KeyError:
			numbers_dict2_new[kk]="XXX"
			chains_dict2_new[kk]="-"
			chains_array2_newkeys.insert(i, kk)
		i +=1
		
	dict2_keys=numbers_dict2_new.keys()
	for kk in chains_array2_newkeys:
		try:
			numbers_dict1_new[kk]
		except KeyError:
			numbers_dict1_new[kk]="XXX"
			chains_dict1_new[kk]="-"
			chains_array1_newkeys.insert(i, kk)

	chains_dict1_keys= chains_dict1_new.keys()
	chains_dict1_keys.sort()
	chains_array1=[] # Array will contain first sequence with deletions
	chains_array2=[] # Array will contain second sequence with deletions


        for kk in chains_array1_newkeys: # At this point all keys are the same in both dictionaries. Now we are filling arrays.
		chains_array1.append(chains_dict1_new[kk])
        for kk in chains_array2_newkeys:		
		chains_array2.append(chains_dict2_new[kk])


	out = []
	pdb_id = os.path.split(pdb_file1)[-1][:-4]
	pdb_id = "%s_%s"%(pdb_id, 1)
	out.append("> %s\n" % (pdb_id))
	
	# Write output in lines 80 characters long, first sequence
	out_new_len = len(chains_array1)
	num_lines= out_new_len/80
	
	for i in range(num_lines+1):
		#out.append("".join([aa for aa in chain_dict1_new[c][80*i:80*(i+1)]]))
		#out.append("\n")
		
		out.append("".join([aa for aa in chains_array1[80*i:80*(i+1)]]))
		out.append("\n")
		
		#out.append("".join([aa for aa in chain_dict1_new[c][80*(i+1):]]))
		#out.append("\n")

	print >>fout, "".join(out)
	out = []
	pdb_id = os.path.split(pdb_file2)[-1][:-4]
	pdb_id = "%s_%s"%(pdb_id, 2)
	out.append("> %s\n" % (pdb_id))

	# Write output in lines 80 characters long, second sequence
	out_new_len = len(chains_array2)
	num_lines= out_new_len/80
	
	for i in range(num_lines+1):
		#out.append("".join([aa for aa in chain_dict2_new[c][80*i:80*(i+1)]]))
		#out.append("\n")
		
		out.append("".join([aa for aa in chains_array2[80*i:80*(i+1)]]))
		out.append("\n")
		
		#out.append("".join([aa for aa in chain_dict2_new[c][80*(i+1):]]))
		#out.append("\n")
	print >>fout, "".join(out)
	fout.close()


# Parse the ATOM and HETATM records in PDB file
# Return dictionaries of atomic names and numbers keyed to serial number of residues
#
def Seq2Atom(pdb):
	res_atomname_dict={} # Set correspondence of Residue serial number and atom names in this residue
	res_atomnumber_dict={} # Set correspondence of Residue serial number and atom serial number 
	atomname_array=[]
	atomnumber_array=[]
	ires_count=0
	atoms = [l for l in pdb if l[0:6] == "ATOM  " or l[0:6] == "HETATM"]
	atomline = atoms[0] # Read first line
	resnum_old = atomline[22:26] # Set first residue number
	if ( atomline[20] == " " ):
		chain_id_old= atomline[21]
	else:
		chain_id_old= atomline[20:22]
	for j in range(len(atoms) ): # Loop over PDB records
		atomline = atoms[j]
		resnum_new = atomline[22:26]
		if ( atomline[20] == " " ):
			chain_id_new= atomline[21]
		else:
			chain_id_new= atomline[20:22]
		if resnum_new != resnum_old or chain_id_new != chain_id_old:
			resnum_old = resnum_new
			chain_id_old = chain_id_new
			ires_count=ires_count+1
			atomname_array=[]
			atomnumber_array=[]
		atomname= atomline[12:16]
		atomnum= j
		atomname_array.append(atomname)# Add atom name
		atomnumber_array.append(atomnum) # Add atom serial number 
		res_atomname_dict[ires_count] = atomname_array
		res_atomnumber_dict[ires_count] = atomnumber_array

	return res_atomname_dict, res_atomnumber_dict
	
# Read Alinment file
# Return dictionary of one letter resudue names keyed to serial number of sequence (sequence 0 and 1)
def ReadAlignFile(alignmentfilename):
	seq={}
	seqline=[]
	iseq = -1
	f = open(alignmentfilename,'r')

        alignfile = f.readlines()
        f.close()
	for i in range (len(alignfile) ):
		strline = alignfile[i]
		if strline[0] != "%" and strline[0] != "#":
			for j in range(len(strline) -1):
				seqline.append(strline[j])
			seq[iseq] = seqline
		else:	
			iseq= iseq +1
			seqline=[]
			
	return seq	

# Find matching atoms in two structures using sequence alingment.
# Also prepare alignment statistics log file that will be used in Pathways site
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

	iseq_len1= 0
	iseq_len2= 0
	icount_first = 0
	icount_second = 0
	iinsertions=0
	idelitions=0
	# Loop over aligned sequences
	for i in range ( len(seq1) ): # Length of seq1 and seq2 is the same at this point (they are aligned)
		atoms_res1=[] # Array of residue atoms from sequence1
		atoms_res2=[] # Array of residue atoms from sequence2
		res1= seq1[i]
		res2= seq2[i]
		control_flag=0
		if res1 != "-" and res2 != "-":
			#print "Case 1 first line",  icount_first, icount_second,
			atoms_res1.append( res_atomname_dict1[icount_first] )
			atoms_res1.append( res_atomnumber_dict1[icount_first] )
			atoms_res2.append( res_atomname_dict2[icount_second] )
			atoms_res2.append( res_atomnumber_dict2[icount_second] )
			icount_first = icount_first+1
			icount_second = icount_second+1
			#print "Case 1",  icount_first, icount_second, res1, res2
			control_flag=1 
		elif res1 != "-" and res2 == "-":
			icount_first = icount_first+1
			iinsertions = iinsertions+1
#			print "Case 2", icount_first, icount_second, res1, res2			
		elif res1 == "-" and res2 != "-":
			icount_second = icount_second+1
			idelitions = idelitions +1
#			print "Case 3",  icount_first, icount_second, res1, res2
		elif res1 == "-" and res2 == "-":
			print "Residues %s - %s are not defined in %s alignment file"%(res1, res2, seqalignfile)
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
			#print icount
		else:
			icount =0
	seq_list=[]
	seq_list.append(seq1)
	seq_list.append(seq2)

	return 	seq_list


def mapAtoms(file_list):
	"""
	
	"""
	extract_seq( file_list[0], file_list[1] )
	seq_list = parseAlignment()
	atomMatch= MatchAtomNumber(file_list, seq_list)
	out=[]

	atom_matchfile ="atommatch.txt"	
	for i in range( len(atomMatch) ):
		out.append(atomMatch[i])
		out.append("\n")
	fout = open(atom_matchfile,'w')
	print >>fout, "".join(out)	
	fout.close()

        print "Finished mapping atoms. Result file:",atom_matchfile
    
	return atom_matchfile
	
if __name__ == "__main__":
    main()
