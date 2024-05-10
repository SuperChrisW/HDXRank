/*
 * ForbidList.cpp
 *
 *  Created on: Jun 21, 2010
 *      Author: dwfarrel
 */

#include "ForbidList.h"
#include "Targeter.h"
#include "ProteinInfo.h"
#include "NeighborTable.h"
#include <algorithm>

using namespace std;

void ForbidList::inserthb( int d, int a ) {
  insert( d, a );
  if ( prot->atom(a).elem() == "O" && (*nt)[a].size() == 1 ) {
    int c = (*nt)[a][0];
    for ( size_t i = 0; i < (*nt)[c].size(); i++ ) {
      int neigh = (*nt)[c][i];
      if ( neigh == a ) continue;
      if ( prot->atom( neigh ).elem() == "O" && (*nt)[neigh].size() == 1 ) {
        insert( d, neigh );
      }
    }
  }
}

void ForbidList::insert( int p1, int p2 ) {
  if ( p2 > p1 ) swap( p1, p2 );
  timestampLookup[pair<int,int>( p1, p2 )] = targeter->getIteration();
}

bool ForbidList::checkforbid( int p1, int p2 ) const {
  if ( p2 > p1 ) swap( p1, p2 );
  //if the pair is found, then we want to forbid the pair. return true.
  return ( timestampLookup.find( pair<int,int>( p1, p2 ) ) != timestampLookup.end() );
}

void ForbidList::tidy() {
  int removethesefromlist = targeter->getIteration() - forbidTimeDuration;
  map<pair<int,int>, int >::iterator it = timestampLookup.begin();
  map<pair<int,int>, int >::iterator temp;
  while ( it != timestampLookup.end() ) {
    if ( it->second <= removethesefromlist ) {
      temp = it;
      temp++;
      timestampLookup.erase( it );
      it = temp;
    }
    else it++;
  }
}
