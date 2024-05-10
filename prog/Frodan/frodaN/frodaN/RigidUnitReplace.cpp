/*
 * RigidUnitReplace.cpp
 *
 *  Created on: Sep 12, 2008
 *      Author: dwfarrel
 */

#include "RigidUnitReplace.h"
#include "RigidUnitSystem.h"
#include "TargetEnergy.h"
#include "mt19937ar.h"
#include "NeighborTable.h"

using namespace std;

//the atommap and cov tables are only needed during the constructor.
//the rigid unit system and target energy are required throughout
//the existence of the RigidUnitReplace object.
RigidUnitReplace::RigidUnitReplace(
    RigidUnitSystem *rigidUnitSystem_,
    const TargetEnergy *targetEnergy_,
    const map<int,int> *atommap,
    const NeighborTable* covTableA,
    const NeighborTable* covTableB ) :
      rigidUnitSystem( rigidUnitSystem_ ),
      targetEnergy( targetEnergy_ )
{

  //make list of replaceable rigid units
  int nRU = rigidUnitSystem->nRigidUnits();
  isReplaceable.assign( nRU, 0 );
  for ( int ru = 0; ru < nRU; ru++ ) {
    if ( checkRigidUnit( ru, atommap, covTableA, covTableB ) ) {
      isReplaceable[ru] = 1;
    }
  }
}

void RigidUnitReplace::initializePQ(
    double rmsdupper,
    double rmsdlower ) {
  double scale = rmsdupper - rmsdlower;
  PQNode pqnode;
  int N = isReplaceable.size();
  for ( int ru = 0; ru < N; ru++ ) {
    if ( !isReplaceable[ru] ) continue;
    pqnode.rmsd = ( scale > 0 ) ? genrand_real2()*scale + rmsdlower : rmsdlower;
    pqnode.ru = ru;
    pq.push( pqnode );
  }
}

bool RigidUnitReplace::checkRigidUnit( int ru,
    const map<int,int> *atommap,
    const NeighborTable* covTableA,
    const NeighborTable* covTableB ) {
  const vector<int>* plistA = &rigidUnitSystem->getPlistFromRU( ru );

  //Make sure each atom in the rigid unit has a target counterpart.
  //Also, look at each covalent bond within the rigid unit, and verify
  //that it is also found in the target covalent bond network
  size_t nP = plistA->size();
  for ( size_t i = 0; i < nP; i++ ) {
    int pA = (*plistA)[i];
    int pB;
    if ( atommap->find( pA ) == atommap->end() ) return false;
    else pB = atommap->find( pA )->second;

    //check the bonds of atom pA
    for ( vector<int>::const_iterator it = (*covTableA)[pA].begin();
          it != (*covTableA)[pA].end(); it++ ) {
      int neighA = *it;
      if ( pA >= neighA ) continue; //avoid double-checking

      //skip bonds that are outside the rigid unit
      if ( find( plistA->begin(), plistA->end(), neighA ) == plistA->end() ) continue;

      //if we arrive here, we have found a covalent bond between atoms in the rigid unit.
      //The question is, does this bond exist in the target structure?
      //if not, then we cannot perform a rigid unit replacement on this rigid unit, so return false.

      int neighB;
      //lookup the neighbor in structure B
      if ( atommap->find( neighA ) == atommap->end() ) {
        return false;
      }
      else {
        neighB = atommap->find( neighA )->second;
      }

      //if the bond does not exist, return false
      if ( find( (*covTableB)[pB].begin(), (*covTableB)[pB].end(), neighB ) ==
             (*covTableB)[pB].end() )
        return false;
    }
  }

  return true;
}

RigidUnitReplace::~RigidUnitReplace() {
}

void RigidUnitReplace::replaceRU( int ru ) {
  if ( !isReplaceable[ru] ) return;

  const vector<int> *ruplist = &rigidUnitSystem->getRUPlistFromRU( ru );
  size_t nRUP = ruplist->size();
  vector<Vec3> newpositions( nRUP );
  for ( size_t i = 0; i < nRUP; i++ ) {
    int p = rigidUnitSystem->getPfromRUP( (*ruplist)[i] );
    bool found = false;
    targetEnergy->getTargetPosition( p, found, newpositions[i] );
    if ( !found ) return;
  }

  fit.setSourceAbsolutePoints( newpositions );
  fit.setTargetAbsolutePoints( rigidUnitSystem->absolutePositions(), *ruplist );
  fit.simpleFit();
  fit.getFitAbsolutePoints( newpositions );

  rigidUnitSystem->setRigidUnit( ru, newpositions );
  isReplaceable[ru] = 0;
}

void RigidUnitReplace::replace( int& count ) {
  count = 0;
  while ( !pq.empty() && targetEnergy->getRMSDtoTarget() < pq.top().rmsd ) {
    replaceRU( pq.top().ru );
    pq.pop();
    count++;
  }
  if ( count ) rigidUnitSystem->update();
}
