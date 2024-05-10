#ifndef NEIGHBORTABLE_H_
#define NEIGHBORTABLE_H_

#include <algorithm>
#include <ostream>
#include <vector>
#include <set>

class NeighborTable
{
public:
  NeighborTable( int Npoints );
  void insert( const NeighborTable& otherNeighborTable );
  void insert( int index1, int index2 );
	virtual ~NeighborTable();

  //shortcut accessors
  //vector<int>& operator[]( size_t i ) { return table1[i]; }
  const std::vector<int>& operator[]( size_t i ) const { return table1[i]; }
  size_t size() const { return table1.size(); }
  const std::vector<int>& getThirdNeighbors_onlyPairsILessThanJ( int i ) const {
    return table3_onlyPairsILessThanJ[i];
  }

  bool isThirdNeighbor( int p1, int p2 ) const {
    //This function is only valid after the user has called the "commit"
    //function once.  That function builds a third-neighbor table.
    //Here, we rely on that table.
    //
    //first, ensure that p1 < p2
    if ( p1 > p2 ) std::swap( p1, p2 );
    std::vector<int>::const_iterator neigh = table3_onlyPairsILessThanJ[p1].begin();
    std::vector<int>::const_iterator end = table3_onlyPairsILessThanJ[p1].end();
    while ( neigh != end && *neigh < p2 ) { neigh++; }
    //when the above while loop is done, we know that either neigh equals end
    //or neigh points to the first element >= p2.
    return neigh != end && *neigh == p2;
  }
  void commit();

  bool isSecondNeighbor( int p1, int p2 ) const {
    //Rely on the fact that the neighbors are listed in sorted order,
    //and that the neighbor table is symmetric (that is, for two
    //points that are neighbors, they are each found in each other's lists)
    //We can look up the neighbors of p1 and look up the neighbors of
    //p2.  If there is at least one common element of both lists,
    //then p1 and p2 are second neighbors.
    if ( p1 == p2 ) return false;
    std::vector<int>::const_iterator n1 = table1[p1].begin();
    std::vector<int>::const_iterator n2 = table1[p2].begin();
    std::vector<int>::const_iterator n1end = table1[p1].end();
    std::vector<int>::const_iterator n2end = table1[p2].end();
    while ( n1 != n1end && n2 != n2end ) {
      if ( *n1 < *n2 ) n1++;
      else if ( *n1 > *n2 ) n2++;
      else return true; // we found a common element, *n1 == *n2
    }
    return false;
  }

  int neighborSearch( int p1, int p2, int maxNeighbor ) const {
    if ( p1 == p2 ) return 0;
    std::set<int> neighborsPrevious;
    std::set<int> neighborsCurrent;
    std::set<int> alreadyVisited;
    std::set<int>::const_iterator neighPrev;
    std::set<int>::const_iterator neighPrevEnd;
    std::vector<int>::const_iterator neighCurrent;
    std::vector<int>::const_iterator neighCurrentEnd;

    neighborsPrevious.insert( p1 );
    alreadyVisited.insert( p1 );
    for ( int L = 1; L <= maxNeighbor; L++ ) {

      //use the current list of neighbors to find the next level of neighbors.
      //we do not add any points that have been already visited,
      //so keep track of visited points
      neighPrev = neighborsPrevious.begin();
      neighPrevEnd = neighborsPrevious.end();
      for ( ; neighPrev != neighPrevEnd; neighPrev++ ) {
        neighCurrent = table1[*neighPrev].begin();
        neighCurrentEnd = table1[*neighPrev].end();
        for ( ; neighCurrent != neighCurrentEnd; neighCurrent++ ) {

          //if the current neighbor equals p2, then we return L, the level
          //at which p1 and p2 are neighbors.
          if ( *neighCurrent == p2 ) return L;

          //only store the current neighbor if it has not already been visited
          if ( alreadyVisited.find( *neighCurrent ) == alreadyVisited.end() ) {
            neighborsCurrent.insert( *neighCurrent );
            alreadyVisited.insert( *neighCurrent );
          }
        }
      }

      swap( neighborsPrevious, neighborsCurrent );
      neighborsCurrent.clear();
    }
    return 0;
  }

  bool isFirstNeighbor( int p1, int p2 ) const {
    //Rely on the fact that the neighbors are listed in sorted order,
    //and that the neighbor table is symmetric (that is, for two
    //points that are neighbors, they are each found in each other's lists)
    std::vector<int>::const_iterator neigh = table1[p1].begin();
    std::vector<int>::const_iterator end = table1[p1].end();
    while ( neigh != end && *neigh < p2 ) { neigh++; }
    //when the above while loop is done, if neigh equals "end",
    //then the search for p2 in the neighbor-list of p1 was unsuccessful.
    //Otherwise, and if *neigh == p2, then the search was successful.
    return neigh != end && *neigh == p2;
  }

  friend std::ostream& operator<< (std::ostream& os, const NeighborTable& neighborTable );
private:
  void setupSecondNeighborTable();
  void setupThirdNeighborTable();
  std::vector< std::set<int> > *table2_;
  std::vector< std::set<int> > *table3_;
  std::vector< std::vector<int> > table1;
  std::vector< std::vector<int> > table3_onlyPairsILessThanJ;
};

#endif /*NEIGHBORTABLE_H_*/