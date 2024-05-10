/*
 * Repulsion.cpp
 *
 *  Created on: Jul 27, 2009
 *      Author: dwfarrel
 */

#include "Repulsion.h"
#include "PairTypeCutoffs.h"
#include "NeighborTable.h"
#include "ProteinInfo.h"
#include <algorithm>
#include <cmath>

using namespace std;

Repulsion::Repulsion(
  const RigidUnitSystem* sys,
  const NeighborTable* nt,
  PairTypeCutoffs* pairTypeCutoffs_ ) :
  testForExclusion( sys, nt ),
  pairTypeCutoffs( pairTypeCutoffs_ )
{
  //build the overriding constraint lookup,
  //starting with all bonded 1-4 pairs inactive ( value -1.0 ).
  //It is not necessary to put the bonded 1-2 and 1-3 pairs, because
  //these are already excluded in the "testForExclusion" function
  int N = sys->nPoints();
  lookupExcludePair_I_LT_J.resize( N );

  //for each atom, make note of 3rd neighbors to exclude
  for ( int i = 0; i < N; i++ ) {
    vector<int>::const_iterator it;
    vector<int>::const_iterator end = nt->getThirdNeighbors_onlyPairsILessThanJ(i).end();
    for ( it = nt->getThirdNeighbors_onlyPairsILessThanJ(i).begin();
          it != end; it++ ) {
      int j = *it;
      lookupExcludePair_I_LT_J[i].insert( j );
    }
  }

  maxCutoff = pairTypeCutoffs->getMaxCutoff();
}

Repulsion::~Repulsion() {
}

void Repulsion::exclude( int i, int j ) {
  if ( i > j ) swap( i, j );
  lookupExcludePair_I_LT_J[i].insert(j);
}

void Repulsion::getCutoff( int p1, int p2, bool& isPairExcluded, double& cutoff ) const {
  //check if self, first, or second neighbors, or same rigid unit.  If so, exclude.
  //This can be done using the TestForExclusion, modified to not exclude 3rd neighbors.
  if ( testForExclusion( p1, p2 ) ) {
    isPairExcluded = true;
    return;
  }

  if ( p1 > p2 ) swap( p1, p2 );
  set<int>::const_iterator it = lookupExcludePair_I_LT_J[p1].find( p2 );
  if ( it != lookupExcludePair_I_LT_J[p1].end() ) {
    isPairExcluded = true;
    return;
  }

  //If not yet excluded, lookup cutoff from pair type table
  isPairExcluded = false;
  cutoff = pairTypeCutoffs->getCutoff( p1, p2 );
}
