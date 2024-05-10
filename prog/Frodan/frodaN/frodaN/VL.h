#ifndef _VL_H_
#define _VL_H_

#include "Vec3.h"
#include "CellList.h"
#include <vector>
#include <map>
#include <cmath>

// class VerletList is used for keeping track of
// pairs of points that are within interacting distance.  It uses
// a smart updating feature that only refreshes its internal lists of
// interacting points when necessary, described below.  The user interacts
// with the VerletList as follows:
//   First, insert points into the VerletList, specifying the interaction
//     cutoff distance of each.
//   Commit the VerletList, which locks the list of points and puts the
//     points into a spanner.  The user here must supply the desired
//     spanner parameters, and also the desired cushion parameter (see below).
//   Update the proximity monitor as points move.
//   Query the VerletList for the set of points that are within the
//     interaction cutoff distance of a particular point.  You will get some
//     extras here, particularly points that were in the "cushion region"
//     (see below).
//
// Here's how the smart updating works.  These details are abstracted away from the user.
// When the user calls the update
// member function, the VerletList queries the Spanner for a list of
// nearby points relative to every point.  But it does not just query out to the
// interaction cutoff of each point.  Instead, it queries out to the cutoff plus
// an extra cushion distance.  For each point, the VerletList stores a list of nearby points,
// which includes any extra points that lie in the cushion region.  As long as points
// do not move more than cushion/2, this stored list of nearby points is guaranteed
// to contain all points within the cutoff distance.
//
// When the user subsequently calls update(), the VerletList will check
// to see how far points have moved.  If points have not moved more than cushion/2,
// then it will not even bother querying the spanner.  These update() calls can be
// very fast.  As points keep moving, at some point one of the points will wander
// more than cutoff/2.  Then, when update() is called, the VerletList knows that
// its lists of nearby points are no longer guaranteed to be accurate, so it re-queries
// the spanner to get up-to-date lists.
//
////////////////////////////////////////////////////////////////////////////////

class VLPoint {
public:
  VLPoint() {}
  VLPoint( int externalID_, const Vec3* position_ ) :
    externalID( externalID_ ),
    position( position_ ),
    savedPosition( *position ) {}
  int externalID;
  const Vec3 *position;
  Vec3 savedPosition;
};

// Requirements of the template classes:
// Cutoff must define the following functions:
//   void getCutoff( int p1, int p2, bool& isPairExcluded, double& cutoff ) const;
//   double getMaxCutoff() const;
//
// and Collect must define
//   void clearPairs();
//   void addPair( int p1, int p2, double cutoff );
template < class Cutoff, class Collect >
class VL
{
public:
  VL( const Cutoff* cutoff, Collect *collect, double trigger_ );
  virtual ~VL();
  void insert( int externalID, const Vec3 *position );
  void update() {
    if ( isTriggerExceeded() ) {
      makeNewList();
    }
  }
  void makeNewList() {
    recordPositions();
    int count;
    findPairs( true, count );
  }
  int countPairs() {
    int count;
    findPairs( false, count );
    return count;
  }
protected:
  const Cutoff* cutoff;
  Collect *collect;
  double trigger;
  CellList cellList;
  std::vector<VLPoint> points;

  void recordPositions();
  bool isTriggerExceeded() const;
  void findPairs( bool doMakeNewList, int& count );
};

template < class Cutoff, class Collect >
VL< Cutoff, Collect >::VL( const Cutoff* cutoff_, Collect *collect_, double trigger_ ) :
  cutoff(cutoff_),
  collect(collect_),
  trigger(trigger_),
  cellList( cutoff->getMaxCutoff() + 2.0*trigger ) {
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   VerletList Destructor.
////////////////////////////////////////////////////////////////////////////////
template < class Cutoff, class Collect >
VL< Cutoff, Collect >::~VL() {
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Insert a point into the VerletList.  Points can be inserted only when
//   setting up the VerletList.  After you call commit(), you cannot
//   insert points anymore.
// Parameters:
//   position - the address of the point's 3-D Vector.  As the point moves,
//     the proximity monitor will always know what its new position is.
//     IMPORTANT - this address must always be valid!  If the 3-D Vector
//     ends up moving to a different location in memory, this will be bad.
//     For example, if the 3-D Vector is stored in a std::vector<>, and then
//     you extend the size of the std::vector, the std::vector may move its
//     contents somewhere else, invalidating all pointers to its contents.
//   cutoffDist - defines the center-to-center distance beyond which points are
//     not neighbors.  When you query the VerletList for this point, it will return points
//     that are within cutoffDist of this point.
//   id - the ID of the point
// Return Value List:
//   none
////////////////////////////////////////////////////////////////////////////////
template < class Cutoff, class Collect >
void VL< Cutoff, Collect >::insert( int externalID, const Vec3* position ) {
  points.push_back( VLPoint( externalID, position ) );
  int internalID = points.size() - 1;
  cellList.insert( internalID, position );
}

template < class Cutoff, class Collect >
void VL< Cutoff, Collect >::recordPositions() {
  const size_t N = points.size();
  for ( size_t i = 0; i < N; i++ ) {
    points[i].savedPosition = *points[i].position;
  }
}

template < class Cutoff, class Collect >
bool VL< Cutoff, Collect >::isTriggerExceeded() const {
  //Compare current point position to the last saved positions.
  //If points have not moved more than trigger distance,
  //then there is no need to update the list
  //of possible contacts
  double trigger2 = trigger*trigger;
  const size_t N = points.size();
  for ( size_t i = 0; i < N; i++ ) {
    double thisMoved2 = points[i].savedPosition.dist2( *points[i].position );
    if ( thisMoved2 > trigger2 ) return true;
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   The spanner is queried for each point, and results stored in the proximity
//   lists.  The query for each point goes out to its cutoff distance, plus the
//   cushion distance.  The isPairDisqualified() function is checked for each
//   pair, to filter out unwanted pairs from the proximity lists.
// Parameters:
//   none
// Return Value List:
//   none
////////////////////////////////////////////////////////////////////////////////
template < class Cutoff, class Collect >
void VL< Cutoff, Collect >::findPairs( bool doMakeNewList, int& count ) {

  if ( doMakeNewList ) collect->clearPairs();
  cellList.update();

  int internalID1 = -1;
  int internalID2 = -1;
  bool valid = false;
  double cushion = 2.0*trigger;
  count = 0;
  do {
    cellList.getNextPair( internalID1, internalID2, valid );
    if ( !valid ) break;

    int externalID1 = points[internalID1].externalID;
    int externalID2 = points[internalID2].externalID;
    if ( externalID1 >= externalID2 ) continue; //avoid self-counting and double-counting

    double cutoffDist;
    bool isPairExcluded;
    cutoff->getCutoff( externalID1, externalID2, isPairExcluded, cutoffDist );
    if ( isPairExcluded ) continue;

    double distToCheck = cutoffDist + cushion;
    Vec3 delta = *points[internalID2].position;
    delta -= *points[internalID1].position;

    if ( std::abs(delta.x) > distToCheck || std::abs(delta.y) > distToCheck || std::abs(delta.z) > distToCheck ) continue;
    if ( delta.norm2() > distToCheck*distToCheck ) continue;

    if ( doMakeNewList ) collect->addPair( externalID1, externalID2, cutoffDist );
    count++;

  } while ( true );
}

#endif
