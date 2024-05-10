#ifndef TESTFOREXCLUSION_H_
#define TESTFOREXCLUSION_H_

#include "RigidUnitSystem.h"
class NeighborTable;
#include <vector>
#include <algorithm>

class TestForExclusion
{
public:
	TestForExclusion( const RigidUnitSystem *rigidUnitSystem_, const NeighborTable *neighborTable );
	virtual ~TestForExclusion() {}
  bool operator()( int p1, int p2 ) const;
private:
  const RigidUnitSystem *rigidUnitSystem;
  std::vector< std::vector< int > > exclusionlist;
};

inline bool TestForExclusion::operator()( int p1, int p2 ) const {
  if ( rigidUnitSystem->doPointsBelongToSameRigidUnit( p1, p2 ) ) return true;
  
  //this function requires that exclusionlist[p1] 
  //and exclusionlist[p2] be sorted.  It will pick the smaller
  //of the two lists and search it.
  return ( exclusionlist[p1].size() < exclusionlist[p2].size() ) ? 
          binary_search( exclusionlist[p1].begin(), exclusionlist[p1].end(), p2 ) :
          binary_search( exclusionlist[p2].begin(), exclusionlist[p2].end(), p1 );
}
#endif /*TESTFOREXCLUSION_H_*/
