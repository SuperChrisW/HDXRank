/*
 * DistConstraint.cpp
 *
 *  Created on: Jul 24, 2009
 *      Author: dwfarrel
 */

#include "DistConstraint.h"
#include <cmath>

using namespace std;

DistConstraint::DistConstraint() {
}

DistConstraint::~DistConstraint() {
}

double DistConstraint::energy() {
  double dist2 = calcDist2();
  if ( doesDistSquaredMeetConstraint( dist2 ) ) return 0;
  double r = sqrt( dist2 ) - cutoff;
  // The energy is 1/2 * k * r^2, where r = dist - cutoff.
  return k*0.5*r*r;
  //An equivalent expression is k*( 0.5*(dist2 + cutoff2) - sqrt(dist2*cutoff2) );
}

void DistConstraint::addToGradient_P(
  vector<Vec3> &dV_dr_P,
  vector<SecondDerivative> &secondDerivative_P ) {

  Vec3 delta = (*positions)[p2];
  delta -= (*positions)[p1];

  double deltax2 = delta.x*delta.x;
  double deltay2 = delta.y*delta.y;
  double deltaz2 = delta.z*delta.z;
  double dist2 = deltax2 + deltay2 + deltaz2;

  if ( doesDistSquaredMeetConstraint( dist2 ) ) return;
  //if ( dist2 < cutoff2 ) return; //for maxdist

  double cutoff_over_dist = cutoff/sqrt(dist2);

  Vec3 dV_dr( delta );
  double factor = k*(cutoff_over_dist-1.0);
  dV_dr *= factor;

  dV_dr_P[p1] += dV_dr;
  dV_dr_P[p2] -= dV_dr;

  SecondDerivative secondDerivative;

  secondDerivative.d2V_dx2 = k*( 1.0 - cutoff_over_dist*( 1.0 - deltax2/dist2 ) );
  secondDerivative.d2V_dy2 = k*( 1.0 - cutoff_over_dist*( 1.0 - deltay2/dist2 ) );
  secondDerivative.d2V_dz2 = k*( 1.0 - cutoff_over_dist*( 1.0 - deltaz2/dist2 ) );
  secondDerivative.d2V_dxdy = k*cutoff_over_dist/dist2*delta.x*delta.y;
  secondDerivative.d2V_dydz = k*cutoff_over_dist/dist2*delta.y*delta.z;
  secondDerivative.d2V_dzdx = k*cutoff_over_dist/dist2*delta.z*delta.x;

  //Below is an approximate form of the second derivative elements,
  //obtained by assuming that the dist is very close to the cutoff,
  //so that cutoff_over_dist is equal to 1.

  //double factor2 = k/dist2;
  //secondDerivative.d2V_dx2 = factor2*deltax2;
  //secondDerivative.d2V_dy2 = factor2*deltay2;
  //secondDerivative.d2V_dz2 = factor2*deltaz2;
  //secondDerivative.d2V_dxdy = factor2*delta.x*delta.y;
  //secondDerivative.d2V_dydz = factor2*delta.y*delta.z;
  //secondDerivative.d2V_dzdx = factor2*delta.z*delta.x;

  secondDerivative_P[p1] += secondDerivative;
  secondDerivative_P[p2] += secondDerivative;
}
