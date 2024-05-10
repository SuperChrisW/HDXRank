/*
 * MinAngleConstraint.cpp
 *
 *  Created on: Jul 20, 2009
 *      Author: dwfarrel
 */

#include "MinAngleConstraint.h"
#include <cmath>
#include "Gradient.h"

using namespace std;

MinAngleConstraint::MinAngleConstraint() {
}

MinAngleConstraint::~MinAngleConstraint() {
}

double MinAngleConstraint::energy() {
  double Eangle = 0.0;
  Vec3 r12;
  Vec3 r32;
  r12 = (*positions)[p1];
  r12 -= (*positions)[p2];
  r32 = (*positions)[p3];
  r32 -= (*positions)[p2];
  double costheta = r12.dot( r32 ) / sqrt( r12.norm2()*r32.norm2() );

  //if angle is within limit, it contributes zero energy
  if ( costheta < maxCosTheta ) return 0;

  double diff = acos(costheta) - minThetaRad;
  Eangle += k/2.0*diff*diff;
  return Eangle;

}

void MinAngleConstraint::addToGradient_P(
    vector<Vec3> &dV_dr_P,
    vector<SecondDerivative> &secondDerivative_P ) {
  Vec3 r12;
  Vec3 r32;
  r12 = (*positions)[p1];
  r12 -= (*positions)[p2];
  r32 = (*positions)[p3];
  r32 -= (*positions)[p2];
  double r12norm2 = r12.norm2();
  double r12norm = sqrt( r12norm2 );
  double r32norm2 = r32.norm2();
  double r32norm = sqrt( r32norm2 );
  double r12dotr32 = r12.dot( r32 );

  double costheta = r12dotr32 / ( r12norm*r32norm );

  //if angle is within limit, it does not contribute
  if ( costheta < maxCosTheta ) return;

  //compute some vector quantities we will need
  Vec3 r12perp_unitvec = r32;
  r12perp_unitvec -= r12 * (r12dotr32/r12norm2);
  r12perp_unitvec /= sqrt( r12perp_unitvec.norm2() );
  Vec3 r32perp_unitvec = r12;
  r32perp_unitvec -= r32 * (r12dotr32/r32norm2);
  r32perp_unitvec /= sqrt( r32perp_unitvec.norm2() );

  //compute dtheta_dr (vectors) for points r1, r2, and r3
  Vec3 dtheta_dr1 = r12perp_unitvec;
  dtheta_dr1 /= -r12norm;
  Vec3 dtheta_dr3 = r32perp_unitvec;
  dtheta_dr3 /= -r32norm;

  Vec3 dtheta_dr2 = dtheta_dr1 + dtheta_dr3;
  dtheta_dr2.x = -dtheta_dr2.x;
  dtheta_dr2.y = -dtheta_dr2.y;
  dtheta_dr2.z = -dtheta_dr2.z;

  //now, for the three points of our angle, r1, r2, and r3,
  //compute the first and second derivatives
  double prefactor = k*(acos(costheta) - minThetaRad);
  dV_dr_P[p1] += dtheta_dr1*prefactor;
  dV_dr_P[p2] += dtheta_dr2*prefactor;
  dV_dr_P[p3] += dtheta_dr3*prefactor;

  SecondDerivative sd1;
  SecondDerivative sd2;
  SecondDerivative sd3;

  sd1.d2V_dx2 = k*dtheta_dr1.x*dtheta_dr1.x;
  sd1.d2V_dy2 = k*dtheta_dr1.y*dtheta_dr1.y;
  sd1.d2V_dz2 = k*dtheta_dr1.z*dtheta_dr1.z;
  sd1.d2V_dxdy = k*dtheta_dr1.x*dtheta_dr1.y;
  sd1.d2V_dydz = k*dtheta_dr1.y*dtheta_dr1.z;
  sd1.d2V_dzdx = k*dtheta_dr1.z*dtheta_dr1.x;

  sd2.d2V_dx2 = k*dtheta_dr2.x*dtheta_dr2.x;
  sd2.d2V_dy2 = k*dtheta_dr2.y*dtheta_dr2.y;
  sd2.d2V_dz2 = k*dtheta_dr2.z*dtheta_dr2.z;
  sd2.d2V_dxdy = k*dtheta_dr2.x*dtheta_dr2.y;
  sd2.d2V_dydz = k*dtheta_dr2.y*dtheta_dr2.z;
  sd2.d2V_dzdx = k*dtheta_dr2.z*dtheta_dr2.x;

  sd3.d2V_dx2 = k*dtheta_dr3.x*dtheta_dr3.x;
  sd3.d2V_dy2 = k*dtheta_dr3.y*dtheta_dr3.y;
  sd3.d2V_dz2 = k*dtheta_dr3.z*dtheta_dr3.z;
  sd3.d2V_dxdy = k*dtheta_dr3.x*dtheta_dr3.y;
  sd3.d2V_dydz = k*dtheta_dr3.y*dtheta_dr3.z;
  sd3.d2V_dzdx = k*dtheta_dr3.z*dtheta_dr3.x;

  secondDerivative_P[p1] += sd1;
  secondDerivative_P[p2] += sd2;
  secondDerivative_P[p3] += sd3;

}
