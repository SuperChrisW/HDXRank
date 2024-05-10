/*
 * ProteinInfo.cpp
 *
 *  Created on: Jul 22, 2009
 *      Author: dwfarrel
 */

#include "ProteinInfo.h"
#include "AmberPrmtop.h"
#include "PDB.h"
#include <cstdlib>

using namespace std;

ProteinInfo::ProteinInfo( const PDB& pdb ) {
  int n = pdb.atomLines.size();
  int last_resSeq = -1;
  char last_chainID = 0;
  char last_insertionCode = 0;
  for ( int i = 0; i < n; i++ ) {
    int this_resSeq = pdb.atomLines[i].resSeq;
    char this_chainID = pdb.atomLines[i].chainID;
    char this_insertionCode = pdb.atomLines[i].iCode;
    if ( this_resSeq != last_resSeq ||
         this_chainID != last_chainID ||
         this_insertionCode != last_insertionCode ) {
      last_resSeq = this_resSeq;
      last_chainID = this_chainID;
      last_insertionCode = this_insertionCode;
      startNewResidue( pdb.atomLines[i].resName );
    }
    string elem = pdb.atomLines[i].element;
    if ( elem.size() == 0 ) {
      cout << "Error extracting element info from PDB file." << endl;
      exit(0);
    }
    addAtomToCurrentResidue( pdb.atomLines[i].strippedName(), elem );
  }
  vector<Resi>(resi_).swap(resi_);
  vector<Atom>(atom_).swap(atom_);
}

ProteinInfo::ProteinInfo( const AmberPrmtop& prmtop ) {
  int n = prmtop.natom;
  int lastresi = -1;
  for ( int i = 0; i < n; i++ ) {
    int thisresi = prmtop.lookupResidueIndexFromAtomIndex[i];
    if ( thisresi != lastresi ) {
      startNewResidue( prmtop.residueLabels[thisresi] );
      lastresi = thisresi;
    }
    string elem;
    elem = prmtop.amberAtomType[i][0]; //this may not always be true...
    addAtomToCurrentResidue( prmtop.atomName[i], elem );
  }
  vector<Resi>(resi_).swap(resi_);
  vector<Atom>(atom_).swap(atom_);
}
