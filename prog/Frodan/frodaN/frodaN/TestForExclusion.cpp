#include "TestForExclusion.h"
#include "RigidUnitSystem.h"
#include "NeighborTable.h"
#include <algorithm>
using namespace std;

TestForExclusion::TestForExclusion( const RigidUnitSystem *rigidUnitSystem_, const NeighborTable *neighborTable ) :
  rigidUnitSystem( rigidUnitSystem_ )
{
  int natoms = rigidUnitSystem->nPoints();

  //this exclusion set is local to this function only.  It has
  //a set, which automatically keeps elements sorted and unique
  vector< set<int> > exclusionset( natoms );

  //For each atom, we wish to build a list of extra atoms to exclude
  //in addition to those in the same rigid unit as the atom.
  //These extra excluded atoms are any covalently bonded first, second, and
  //third neighbors that are in a different rigid units from the atom.

  //Initially, the exclude list gets all first and second neighbors.
  //But then we go through and remove from the list any that are in
  //the same rigid unit (as this exclude list is only meant to
  //contain the "extra" atoms to exclude)
  for ( int atom = 0; atom < natoms; atom++ ) {
    set<int> *tempset = &exclusionset[atom];
    const vector<int> *neighlist1;
    const vector<int> *neighlist2;
    //const vector<int> *neighlist3;
    int neigh1;
    int neigh2;
    //int neigh3;
    size_t i;
    size_t j;
    //size_t k;
    neighlist1 = &(*neighborTable)[atom];
    for ( i = 0; i < neighlist1->size(); i++ ) {
      neigh1 = (*neighlist1)[i];
      tempset->insert( neigh1 );
      neighlist2 = &(*neighborTable)[neigh1];
      for ( j = 0; j < neighlist2->size(); j++ ) {
        neigh2 = (*neighlist2)[j];
        tempset->insert( neigh2 );
        /*
        neighlist3 = &(*neighborTable)[neigh2];
        for ( k = 0; k < neighlist3->size(); k++) {
          neigh3 = (*neighlist3)[k];
          tempset->insert( neigh3 );
        }
        */
      }
    }

    set<int>::iterator it = tempset->begin();
    set<int>::iterator end = tempset->end();
    set<int>::iterator it_remove;
    while ( it != end ) {
      if ( rigidUnitSystem->doPointsBelongToSameRigidUnit( atom, *it ) ) {
        it_remove = it++;
        tempset->erase( it_remove );
      }
      else it++;
    }

  }

  //copy these sets over to the permanent storage,
  //converting the sets to vectors
  exclusionlist.resize( natoms );
  for ( int atom = 0; atom < natoms; atom++ ) {
    exclusionlist[atom].assign( exclusionset[atom].begin(), exclusionset[atom].end() );
  }

}
