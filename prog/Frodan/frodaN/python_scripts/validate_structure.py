#!/usr/bin/env python

# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
validateStructure.py
Checks if PDB file has Alternative Atoms Positions and TER records. 
If PDB has Alternative Atoms Positions script will output only one position of atom with atom's highest occupancy.
Script also  extracts only first model in multi-model PDB file.
"""
__author__ = "Kirill Speranskiy"
__date__ = "05-12-09"

import os
import sys


def checkChains(atoms): # Check if PDB file has TER records to separate the chains
	chain_id_old = ""
	atomline_old = ""
	atoms_checked = []
	ter_record = ""
	resnum_old = 0
	for j in range( len(atoms) ): # Loop over PDB records
		atomline = atoms[j]
		if (atomline[0:6] == "ATOM  " or atomline[0:6] == "HETATM"):
			chain_id = atomline[20:22]
			chain_id_old = atomline_old[20:22]
			res_name = atomline[17:20]
			res_name_old = atomline_old[17:20]
			resnum = int( atomline[23:26] )
			if (chain_id != chain_id_old):
				before_atomline= atoms[j-1]
				if(before_atomline[0:3] != "TER"  ):
					if ( before_atomline[0:6] == "ATOM  " or before_atomline[0:6] == "HETATM" ):
						ter_record = "TER\n"
				chain_id_old = chain_id
			if (chain_id == "  " and chain_id_old == "  "):
				if (res_name != res_name_old and resnum <= resnum_old):
					before_atomline= atoms[j-1]
					if(before_atomline[0:3] != "TER"):
						ter_record = "TER\n"
			atomline_old = atomline
			resnum_old = resnum
		if (ter_record != ""): # Add TER record
			atoms_checked.append(ter_record)
			ter_record = ""
		atoms_checked.append(atomline)
		
	return atoms_checked
		
def checkWaterAndHydrogen(atomline):
        result=0
	try:
		if (atomline[17:21] != "TIP3" and atomline[17:21] != "TIP4" and atomline[17:21] != "OPLC"):
			atom_name1 = atomline[12]
			atom_name2 = atomline[13]
			if(atom_name1 != "H" and atom_name2 != "H"):
				result =1
				return result
			else:
				result =0
				return result
			
	except:	
		return result
	
def checkAtoms(res_array):
	atom_dict={}
	occupancy1= 0
	occupancy2= 0
	occupancy_letter_new = ""
	occupancy_letter_old = ""	
	for i in range ( len(res_array) ):
		atom = res_array[i]
		atom_name1 = atom[12:16]
		atom_saved = atom_dict.get( atom_name1 )
		
		if ( atom_saved ):
			atom_name2 = atom_saved[12:16]
			if (atom_name1 == atom_name2):
				try:
					occupancy1= float( atom[54:60] )
					occupancy2= float( atom_saved[54:60] )
				except:
					occupancy1= 0
					occupancy2= 0
				if (occupancy1 > occupancy2):						
					atom_dict[atom_name1] = atom
		else:
			atom_dict[atom_name1] = atom

	# Check if all atoms have the same alternate location indicator.
	# If they do not than create new atom_dict and take the first model (e.g., atom variant A)
	alter_pos_counter= 0 
	for index, aa in enumerate(atom_dict):
		atom = atom_dict[aa]
		occupancy_letter_new = atom[16]
		if ( occupancy_letter_new != occupancy_letter_old):
		     alter_pos_counter +=1
		     occupancy_letter_old = occupancy_letter_new
		if (occupancy_letter_new == " "): # Do not check if residue has no alternative positions
			break
	if (alter_pos_counter >1):
		atom_dict={} # take the first model if the alternate location indicator was not the same
		for i in range ( len(res_array) ):
			atom = res_array[i]
			atom_name1 = atom[12:16]
			atom_saved = atom_dict.get( atom_name1 )		
			if ( atom_saved == None):
				atom_dict[atom_name1] = atom	
			
	# Create new array of atoms
	res_arr_new=[]
	for i in range ( len(res_array) ):
		atom = res_array[i]
		atom_name1 = atom[12:16]
		if( len(atom_dict[atom_name1]) > 0 ):
			res_arr_new.append(atom_dict[atom_name1])
			atom_dict[atom_name1]=""

	return res_arr_new
	

def validate(pdb_data, hetatm):
	"""
	Check if PDB structure has proper formating.  
	"""
       
	conect_records = [l for l in pdb_data if l[0:6] == "CONECT"]
	end_record = [l for l in pdb_data if l[0:6] == "END   "]
	pdb_data_init=[]

	if(hetatm):	
		pdb_data_init = [l for l in pdb_data if  l[0:6] != "CONECT" and l[0:6] != "END   "]
	else:
		pdb_data_init = [l for l in pdb_data if l[0:6] != "HETATM" and l[0:6] != "CONECT" and l[0:6] != "END   " ]
	atoms_all = checkChains(pdb_data_init)
	outArray=[]
	alt_location =""
	ind=0
	occupancy_min =0
	resID_old=0
	res_array=[] # Array of atoms from one residue
	for j in range( len(atoms_all) ): # Loop over PDB records
		atomline = atoms_all[j]

		if(atomline[0:6] == "ATOM  " or atomline[0:6] == "HETATM" ):
			#if( checkWaterAndHydrogen(atomline) ):

				atom_name1 = atomline[12:16]
				resID1= atomline[21:26]
				
				if ( resID1 == resID_old ):
					res_array.append(atomline)
										
				else:
					res_array_checked = checkAtoms(res_array)
					for k in range (len (res_array_checked ) ):
						outArray.append(res_array_checked[k])

      					res_array= []	
					res_array.append(atomline)
					resID_old= resID1
				if ( j ==  (len(atoms_all)-1) ):
					if ( len(res_array) > 0 ):
						res_array_checked = checkAtoms(res_array)
						for k in range (len (res_array_checked ) ):
							outArray.append(res_array_checked[k])
						res_array= []
		elif(  atomline[0:3] == "TER" or j == len(atoms_all)-1 ): # For the last line in the file
			if ( len(res_array) > 0 ):
				res_array_checked = checkAtoms(res_array)
				for k in range (len (res_array_checked ) ):
					outArray.append(res_array_checked[k])
			       	res_array= []

			outArray.append(atomline)

		# Extract only first model in multi-model PDB file
		try:
			if (atomline[0:6] == "ENDMDL"):
				break
		except:
			print atomline + " is TOO short"
			
	if ( len(conect_records) != 0):
		for i in range( len(conect_records) ):
			outArray.append(conect_records[i])
	if ( len(end_record) != 0):	
		outArray.append(end_record[0])

        print "validateStructure finished"
	return outArray



