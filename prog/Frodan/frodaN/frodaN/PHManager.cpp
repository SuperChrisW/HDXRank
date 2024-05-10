/*
 * PHManager.cpp
 *
 *  Created on: Aug 10, 2009
 *      Author: dwfarrel
 */


#include "PHManager.h"
#include "ProteinInfo.h"
#include "ForbidList.h"
#include <cmath>
#include <algorithm>
#include <string>
#include <set>

using namespace std;

VerletCutoffsHydrophobicPairs::VerletCutoffsHydrophobicPairs( const ProteinInfo& prot, const NeighborTable* nt_, const vector<char>* isNuc_ ) :
  nt( nt_ ),
  isNuc( isNuc_ ) {
  int natoms = prot.natoms();
  lookupResIndexFromAtom.resize( natoms );
  for ( int i = 0; i < natoms; i++ ) {
    lookupResIndexFromAtom[i] = prot.atom(i).resi().index();
  }
}

VerletCutoffsHydrophobicPairs::~VerletCutoffsHydrophobicPairs() {}

void VerletCutoffsHydrophobicPairs::getCutoff( int p1, int p2, bool& isPairExcluded, double& cutoff ) const {
  cutoff = 5.0;

  //Logic added to handle nucleic acids following
  //work of Holger Gohlke and Simone Fuller
  if ( (*isNuc)[p1] && (*isNuc)[p2] ) {
    isPairExcluded =
      nt->isFirstNeighbor( p1, p2) ||
      nt->isSecondNeighbor( p1, p2 ) ||
      nt->isThirdNeighbor( p1, p2 );
  }
  else {
    isPairExcluded =
      nt->isFirstNeighbor( p1, p2) ||
      nt->isSecondNeighbor( p1, p2 ) ||
      nt->isThirdNeighbor( p1, p2 ) ||
      lookupResIndexFromAtom[p1] == lookupResIndexFromAtom[p2];
  }
}

double VerletCutoffsHydrophobicPairs::getMaxCutoff() const { return 5.0; }



PHManager::PHManager(
  const ProteinInfo& prot_,
  const NeighborTable& nt_,
  const vector<Vec3>& coords_,
  PHContainer& phobes_ ) :
    prot( prot_ ),
    nt(nt_),
    coords(coords_),
    phobes( phobes_ ),
    Lnuc( 3.55 ),
    k( 10.0 ),
    L( 3.9 ),
    nAdded_(0),
    nRemoved_(0),
    forbidlist( NULL ),
    carefulChecking( true )
{

  set<string> hydrophobicSet;
  hydrophobicSet.insert("LEU");
  hydrophobicSet.insert("ILE");
  hydrophobicSet.insert("VAL");
  hydrophobicSet.insert("PHE");
  hydrophobicSet.insert("TRP");
  hydrophobicSet.insert("MET");
  hydrophobicSet.insert("ALA");
  hydrophobicSet.insert("TYR");

  set<string> proteinSet;
  proteinSet.insert( "LEU" );
  proteinSet.insert( "ALA" );
  proteinSet.insert( "VAL" );
  proteinSet.insert( "MET" );
  proteinSet.insert( "ILE" );
  proteinSet.insert( "TRP" );
  proteinSet.insert( "PHE" );
  proteinSet.insert( "SER" );
  proteinSet.insert( "THR" );
  proteinSet.insert( "CYS" );
  proteinSet.insert( "ASN" );
  proteinSet.insert( "GLN" );
  proteinSet.insert( "TYR" );
  proteinSet.insert( "HIS" );
  proteinSet.insert( "ASP" );
  proteinSet.insert( "GLU" );
  proteinSet.insert( "LYS" );
  proteinSet.insert( "ARG" );
  proteinSet.insert( "PRO" );
  proteinSet.insert( "GLY" );


  //to detect nucleic acids, a residue must contain all of the
  //following ribose atoms
  set<string> requiredRiboseAtoms;
  requiredRiboseAtoms.insert("C5\'");
  requiredRiboseAtoms.insert("C4\'");
  requiredRiboseAtoms.insert("C3\'");
  requiredRiboseAtoms.insert("C2\'");
  requiredRiboseAtoms.insert("C1\'");
  requiredRiboseAtoms.insert("O4\'");

  //here is the same set as above, with old atom names
  set<string> requiredRiboseAtomsOld;
  requiredRiboseAtomsOld.insert("C5*");
  requiredRiboseAtomsOld.insert("C4*");
  requiredRiboseAtomsOld.insert("C3*");
  requiredRiboseAtomsOld.insert("C2*");
  requiredRiboseAtomsOld.insert("C1*");
  requiredRiboseAtomsOld.insert("O4*");

  //to detect nucleic acids, a residue must contain all of the
  //following base atoms (it will contain more than these, of course,
  //but this is the minimal required amount)
  set<string> requiredBaseAtoms;
  requiredBaseAtoms.insert("C2");
  requiredBaseAtoms.insert("C4");
  requiredBaseAtoms.insert("C6");

  //carbons in BASES of nucleic acids are treated differently, so here is the list
  //of base carbons.
  set<string> nucBaseNames;
  nucBaseNames.insert("C2");
  nucBaseNames.insert("C4");
  nucBaseNames.insert("C5");
  nucBaseNames.insert("C6");
  nucBaseNames.insert("C8");

  //prepare the verlet list
  verletCutoffs = new VerletCutoffsHydrophobicPairs( prot, &nt, &isNucAtom );
  verletCollector = new VerletCollectorHydrophobicPairs;
  verlet = new VerletList( verletCutoffs, verletCollector, 1.0 );

  int natoms = coords.size();
  int nresi = prot.nresi();

  //initialize isNucAtom to all zeros (false)
  isNucAtom.resize( natoms, 0 );

  //determine whether each residue is a nucleic acid or not.
  //If so, mark all the atoms as nucleic acid atoms.
  for ( int r = 0; r < nresi; r++ ) {

    //load all the atom names of this residue into a set
    Resi::const_iterator resbegin = prot.resi(r).begin();
    Resi::const_iterator resend = prot.resi(r).end();
    set<string> atomnames;
    atomnames.clear();
    for ( Resi::const_iterator atom = resbegin; atom != resend; atom++ ) {
      atomnames.insert( atom->name() );
    }

    //the residue is a nucleic acid if all the names in nucleotideSet are
    //found in res
    bool isNucResi =
      includes( atomnames.begin(), atomnames.end(), requiredBaseAtoms.begin(), requiredBaseAtoms.end() ) &&
      ( includes( atomnames.begin(), atomnames.end(), requiredRiboseAtoms.begin(), requiredRiboseAtoms.end() ) ||
        includes( atomnames.begin(), atomnames.end(), requiredRiboseAtomsOld.begin(), requiredRiboseAtomsOld.end() ) );


    //if it is a nucleic acid, mark all the atoms as nucleic acid atoms
    if ( isNucResi ) {
      for ( Resi::const_iterator atom = resbegin; atom != resend; atom++ ) {
        isNucAtom[atom->index()] = 1;
      }

    }
  }

  isNucBaseAtom.resize( natoms, 0 );
  isAtomInHydrophobicVerletList.resize( natoms, 0 );
  for ( int i = 0; i < natoms; i++ ) {
    string elem = prot.atom(i).elem();
    string name = prot.atom(i).name();
    string resi = prot.atom(i).resi().name();

    //proteins
    if ( proteinSet.find( resi ) != proteinSet.end() ) {
      if ( hydrophobicSet.find( resi ) != hydrophobicSet.end() ) {
        if ( ( elem == "C" || elem == "S" ) &&
             name != "C" &&
             name != "CA"  ) {
             // !( ( resi == "PHE" || resi == "TYR" ) && ( name == "CZ" || name == "CD1" || name == "CD2" ) ) )
          verlet->insert( i, &coords[i] );
          isAtomInHydrophobicVerletList[i] = 1;
	}
      }
    }

    //nucleotides
    else if ( isNucAtom[i] ) {
      if ( elem == "C" ) {
        verlet->insert( i, &coords[i] );
        isAtomInHydrophobicVerletList[i] = 1;
        isNucBaseAtom[i] = nucBaseNames.find( name ) != nucBaseNames.end();
      }
    }

    //other groups
    else if ( ( elem == "C" || elem == "S" ) && hasOnlyCHneighbors( i ) ) {
      verlet->insert( i, &coords[i] );
      isAtomInHydrophobicVerletList[i] = 1;
    }
  }
  verlet->makeNewList();
}

  PHManager::~PHManager() {
  delete verletCutoffs;
  delete verletCollector;
  delete verlet;
}

bool PHManager::hasOnlyCHneighbors( int atom ) {
  for ( size_t i = 0; i < nt[atom].size(); i++ ) {
    int neigh = nt[atom][i];
    if ( prot.atom(neigh).elem() != "C" &&
         prot.atom(neigh).elem() != "H" ) return false;
  }
  return true;
}

bool PHManager::isPairExcluded( int p1, int p2 ) const {
  int N = isAtomInHydrophobicVerletList.size();
  if ( 
    !( p1 < N && p2 < N &&
       isAtomInHydrophobicVerletList[p1] && 
       isAtomInHydrophobicVerletList[p2] ) 
  ) return true;

  bool isPairExcluded;
  double cutoff;
  verletCutoffs->getCutoff( p1, p2, isPairExcluded, cutoff );
  return isPairExcluded;
}

void PHManager::tighten() {
  for ( PHContainer::iterator ph = phobes.begin(); ph != phobes.end(); ph++ ) {
    tightenConstraint( *ph );
  }
}

void PHManager::findnew( bool breakable ) {
  verlet->update();
  VerletCollectorHydrophobicPairs::const_iterator it;
  for ( it = verletCollector->begin(); it != verletCollector->end(); it++ ) {
    bool success;
    addConstraint( it->p1, it->p2, success, breakable );
  }
}

void PHManager::addConstraint( int p1, int p2, bool& success, bool breakable ) {
  success = false;

  if ( forbidlist && forbidlist->checkforbid( p1, p2 ) ) return;

  //if this constraint already exists, do not add it again.
  if ( phobes.find( p1, p2 ) != phobes.end() ) { success = false; return; }

  if ( carefulChecking ) {
    double dist = sqrt( ( coords[p1] - coords[p2] ).norm2() );

    //special handling for nucleic acids
    if ( isNucAtom[p1] && isNucAtom[p2] ) {
      if ( dist > Lnuc ) { success = false; return; }
      // if both are bases, check to see if a constraint already exists for the base pair.
      // if so, then if the new constraint is shorter remove the old constraint and change
      //   the entry for the base pair to reflect the new constraint.
      // If there is no constraint for this base pair, assign it to the base pair.
      if ( isNucBaseAtom[p1] && isNucBaseAtom[p2] ) {
        NucBasePair bp;
        bp.base[0] = prot.atom(p1).resi().index();
        bp.base[1] = prot.atom(p2).resi().index();
        map<NucBasePair,HydrophobicPair>::iterator it = mapBasePairToPhobicPair.find( bp );
        if ( it != mapBasePairToPhobicPair.end() ) {
          int otherp1 = it->second.p1;
          int otherp2 = it->second.p2;
          double distother = sqrt( ( coords[otherp1] - coords[otherp2] ).norm2() );
          if ( dist < distother ) {
            phobes.erase( otherp1, otherp2 );
            it->second.p1 = p1;
            it->second.p2 = p2;
          }
          else { success = false; return; }
        }
        else {
          mapBasePairToPhobicPair[bp] = HydrophobicPair( p1, p2 );
        }
      }
    }

    //regular handling (not nucleic acids)
    //Only add the constraint if the distance is closer than L.
    //Note that the actual constraint distance will be initialized to dist + 0.5.
    else if ( dist > L ) { success = false; return; }
  }

  //all cases (including nucleic acids)
  PHConstraint hydrophobicConstraint;
  nAdded_++;
  initializeConstraint( p1, p2, hydrophobicConstraint );
  if ( !breakable ) hydrophobicConstraint.makeUnbreakable();
  phobes.insert( p1, p2, hydrophobicConstraint );
  success = true;
}

void PHManager::initializeConstraint( int p1, int p2, PHConstraint& hydrophobicConstraint ) {
  //initialize atoms
  hydrophobicConstraint.setPoints( &coords, p1, p2 );

  //initialize dist
  double dist = sqrt( ( coords[p1] - coords[p2] ).norm2() );
  hydrophobicConstraint.setCutoff( dist + 0.5 );

  //initialize k
  hydrophobicConstraint.setk( k  );
}

void PHManager::tightenConstraint( PHConstraint& hydrophobicConstraint ) {
  if ( !hydrophobicConstraint.isBreakable() ) return;

  //tighten the distance constraint
  double constraintDist = hydrophobicConstraint.getCutoff();
  double newConstraintDist = hydrophobicConstraint.calcDist() + 0.5;
  if ( newConstraintDist < constraintDist ) {
    hydrophobicConstraint.setCutoff( newConstraintDist );
  }
}

void PHManager::breakAllBreakable() {

  PHContainer::iterator it = phobes.begin();
  while ( it != phobes.end() ) {
    if ( it->isBreakable() ) {

      //special handling for nucleic acids
      //remove constraint entry from base pair lookup, if it's there
      int p1 = it->getp1();
      int p2 = it->getp2();
      if ( isNucBaseAtom[p1] && isNucBaseAtom[p2] ) {
        NucBasePair bp;
        bp.base[0] = prot.atom(p1).resi().index();
        bp.base[1] = prot.atom(p2).resi().index();
        mapBasePairToPhobicPair.erase( bp );
      }

      //remove constraint
      phobes.erase( it );
      nRemoved_++;
      //DO NOT advance the iterator "it".  The erase operation
      //of ConstraintContainer replaces the current constraint
      //with the final constraint.  So, the "next" constraint
      //to be checked is the one that "it" is currently pointing to.
    }
    else it++;
  }

}
