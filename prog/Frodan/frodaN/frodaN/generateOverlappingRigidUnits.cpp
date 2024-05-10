#include "generateOverlappingRigidUnits.h"
#include "NeighborTable.h"
#include <map>
#include <set>
using namespace std;

void generateOverlappingRigidUnits( 
    const vector<int>& mapPtoRC,
    const NeighborTable& neighborTable,
    std::vector< std::vector<int> >& mapRUtoP_vector)
{
  int nP = mapPtoRC.size();

  map< int, set< int > > mapRCtoP;
  for ( int p = 0; p < nP; p++ ) {
    int rc = mapPtoRC[p];
    mapRCtoP[rc].insert( p );
  }

  // Now extend each rigid cluster to include bonded first neighbors.
  // Iterate over each rigid cluster...
  for ( map< int, set< int > >::iterator it = mapRCtoP.begin();
        it != mapRCtoP.end(); )
  {
    set<int> *pointset = &it->second;
    //if this rigid cluster has only one point,
    //and if that one point has exactly one neighbor,
    //then we remove the rigid cluster.
    //We do this because otherwise after extending
    //all rigid units by one neighbor, we would end up 
    //with a redundant two-point rigid unit
    map< int, set< int > >::iterator next_it;    
    if ( pointset->size() == 1 ) {
      int p = *pointset->begin();
      if ( neighborTable[p].size() == 1 ) {
        next_it = it;
        next_it++;
        mapRCtoP.erase( it );
        it = next_it;
        continue;
      }
    }

    //Identify atoms that should be added to this rigid unit.
    //These are any first neighbors of any of the points
    //in the rigid unit.
    //Collect the atoms to be added in the variable atomsToAdd.
    set<int> atomsToAdd;
    atomsToAdd.clear();
    set<int>::iterator pointiter;
    for ( pointiter = pointset->begin();
          pointiter != pointset->end();
          pointiter++ ) {
      int p = *pointiter;
      
      //collect all of p's first neighbors
      const vector<int> *neighborsOfP = &neighborTable[p];
      atomsToAdd.insert( neighborsOfP->begin(), neighborsOfP->end() );
    }

    //Now add the new atoms to the rigid unit
    for ( pointiter = atomsToAdd.begin();
          pointiter != atomsToAdd.end();
          pointiter++ )
    {
      pointset->insert( *pointiter );
    }
    
    it++;
  }

  //clean up.
  mapRUtoP_vector.resize( mapRCtoP.size() );
  int ru = 0;
  for ( map< int, set< int > >::iterator it = mapRCtoP.begin();
        it != mapRCtoP.end();
        it++ )
  {
    const set<int> *plist = &it->second;
    mapRUtoP_vector[ru].resize( plist->size() );
    copy( plist->begin(), plist->end(), mapRUtoP_vector[ru].begin() );
    ru++;
  }
  
}
