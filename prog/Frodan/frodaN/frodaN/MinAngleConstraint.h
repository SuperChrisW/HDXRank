/*
 * MinAngleConstraint.h
 *
 *  Created on: Jul 20, 2009
 *      Author: dwfarrel
 */

#ifndef MINANGLECONSTRAINT_H_
#define MINANGLECONSTRAINT_H_

#include "Vec3.h"
#include <vector>
#include <cmath>
#include "Gradient.h"
#include "Energy.h"

class MinAngleConstraint :
  public EnergyTerm,
  public GradientTerm_P {
public:
  MinAngleConstraint();
  MinAngleConstraint( const std::vector<Vec3>* positions_, double k_, int p1_, int p2_, int p3_, double minThetaRad_ ) :
    positions(positions_),
    k(k_),
    p1(p1_),
    p2(p2_),
    p3(p3_),
    minThetaRad(minThetaRad_),
    maxCosTheta(cos(minThetaRad_)) {}
  virtual ~MinAngleConstraint();

  double calcCosTheta() const {
    Vec3 r12;
    Vec3 r32;
    r12 = (*positions)[p1];
    r12 -= (*positions)[p2];
    r32 = (*positions)[p3];
    r32 -= (*positions)[p2];
    return r12.dot( r32 ) / sqrt( r12.norm2()*r32.norm2() );
  }
  double calcAngleStretch() const {
    return acos(calcCosTheta()) - minThetaRad;
  }
  double energy();
  void addToGradient_P(
    std::vector<Vec3> &dV_dr_P,
    std::vector<SecondDerivative> &secondDerivative_P );

  double getk() const { return k; }
  int getp1() const { return p1; }
  int getp2() const { return p2; }
  int getp3() const { return p3; }
  double getMinAngleRad() const { return minThetaRad; }
  void setk( double k_ ) { k = k_; }
  void setPoints( const std::vector<Vec3>* positions_, int p1_, int p2_, int p3_ ) {
    positions = positions_; p1 = p1_; p2 = p2_; p3 = p3_; }
  void setMinAngleRad( double minThetaRad_ ) {
    minThetaRad = minThetaRad_;
    maxCosTheta = cos(minThetaRad);
  }

protected:
  const std::vector<Vec3>* positions;
  double k;
  int p1;
  int p2;
  int p3;
  double minThetaRad;
  double maxCosTheta;
};

#endif /* MINANGLECONSTRAINT_H_ */
