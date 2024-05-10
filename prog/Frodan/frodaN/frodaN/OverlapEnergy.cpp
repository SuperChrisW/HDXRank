/*
 * OverlapEnergy.cpp
 *
 *  Created on: Jul 24, 2009
 *      Author: dwfarrel
 */

#include "OverlapEnergy.h"
#include "RigidUnitSystem.h"
#include <cmath>
#include <limits>
#include <iostream>

using namespace std;

OverlapEnergy::OverlapEnergy( RigidUnitSystem* sys_, Repulsion* repulsion ) :
    sys( sys_ ),
    verlet( repulsion, this, 1.0 ),
    k( 10.0 ) {

  positions = &sys->meanPositions();
  sys->registerObserver( this );

  //setup verlet list
  const size_t nPoints = sys->nPoints();
  for ( size_t p=0; p<nPoints; p++ ) {
    verlet.insert( p, &sys->meanPositions(p) );
  }
  //don't actually get the list of pairs yet.
  //First, we will reserve space.
  //Ask verlet list for size of list.  Reserve 120% of that size.
  //if ever the size climbs such that the capacity is only 110% of the size,
  //increase reserve to 120% of the new size.
  int Npairs = verlet.countPairs();
  size_t capacity = static_cast<size_t>( 1.2*Npairs );
  constraints.reserve( capacity );
  sizeTrigger = static_cast<size_t>( 1.1*Npairs );

  verlet.makeNewList();
}

OverlapEnergy::~OverlapEnergy() {
}

double OverlapEnergy::energy() {
  double E = 0;
  const size_t n = constraints.size();
  for ( size_t i = 0; i < n; i++ ) {
    E += constraints[i].energy();
  }
  return E;
}

void OverlapEnergy::addToGradient_P(
    std::vector<Vec3> &dV_dr_Point,
    std::vector<SecondDerivative> &secondDerivative_Point ) {
  const size_t n = constraints.size();
  for ( size_t i = 0; i < n; i++ ) {
    constraints[i].addToGradient_P( dV_dr_Point, secondDerivative_Point );
  }
}

double OverlapEnergy::mismatch() {
  double maxOverlap = 0.0;
  if ( verbose ) {
    cout << "Overlaps:" << endl;
    cout << " point1 point2 overlapDist" << endl;
  }
  int i_worst = 0;
  const size_t n = constraints.size();
  for ( size_t i = 0; i < n; i++ ) {
    double overlapDist = constraints[i].getCutoff() - constraints[i].calcDist();
    if ( overlapDist > maxOverlap ) {
      maxOverlap = overlapDist;
      i_worst = i;
    }
    if ( verbose && overlapDist > 0.1 ) {
      cout << "  " << constraints[i].getp1() << " " << constraints[i].getp2() << " " << overlapDist << endl;
    }
  }
  if ( verbose && maxOverlap > numeric_limits<double>::epsilon() ) {
    double sep = constraints[i_worst].calcDist();
    cout << "Max Overlap: point1 " << constraints[i_worst].getp1() << " point2 " << constraints[i_worst].getp2() <<
    " separation " << sep <<
    " overlapDist " << maxOverlap <<
    " constraint " << constraints[i_worst].getCutoff() << endl;
  }
  return maxOverlap;
}
