/*
 * HBConstraint.h
 *
 *  Created on: Jul 20, 2009
 *      Author: dwfarrel
 */

#ifndef HBCONSTRAINT_H_
#define HBCONSTRAINT_H_

#include <cmath>
#include "DistConstraint.h"
#include "MinAngleConstraint.h"
#include <vector>
#include "Vec3.h"
#include "Energy.h"
#include "Gradient.h"
#include <iostream>

using namespace std;

class HBConstraintAbstract :
  public EnergyTerm,
  public GradientTerm_P {
public:
  HBConstraintAbstract();
  virtual ~HBConstraintAbstract();

  double calcDistHA() const { return distHA.calcDist(); }
  double calcDistSquaredHA() const { return distHA.calcDist2(); }
  double calcAngleRadDHA() const { return acos(angleDHA.calcCosTheta()); }

  virtual double energy() = 0;
  virtual void addToGradient_P(
    std::vector<Vec3> &dV_dr_P,
    std::vector<SecondDerivative> &secondDerivative_P ) = 0;

  void setMaxDistHA( double dist ) { distHA.setCutoff( dist ); }
  void setMinAngleRadDHA( double theta ) { angleDHA.setMinAngleRad( theta ); }
  void setkDistHA( double k ) { distHA.setk( k ); }
  void setkAngleDHA( double k ) { angleDHA.setk( k ); }
  void makeUnbreakable() { isBreakable_ = false; }
  int getd () const { return angleDHA.getp1(); }
  int geth () const { return angleDHA.getp2(); }
  int geta () const { return angleDHA.getp3(); }
  bool isBreakable() const { return isBreakable_; }
  double getConstraintMaxDistHA() const { return distHA.getCutoff(); }
  double getConstraintMinAngleRadDHA() const { return angleDHA.getMinAngleRad(); }
protected:
  bool isBreakable_;
  MaxDistConstraint distHA;
  MinAngleConstraint angleDHA;
};

class HBConstraint : public HBConstraintAbstract {
public:
  HBConstraint();
  virtual ~HBConstraint();
  HBConstraint( const std::vector<Vec3>* positions, int d, int h, int a );
  void setAtoms( const std::vector<Vec3>* positions, int d, int h, int a ) {
    distHA.setPoints( positions, h, a );
    angleDHA.setPoints( positions, d, h, a );
  }
  virtual double energy();
  virtual void addToGradient_P(
    std::vector<Vec3> &dV_dr_P,
    std::vector<SecondDerivative> &secondDerivative_P );
};

class BBHBConstraint : public HBConstraintAbstract {
public:
  BBHBConstraint();
  virtual ~BBHBConstraint();
  BBHBConstraint( const std::vector<Vec3>* positions, int d, int h, int a, int b );
  void setAtoms( const std::vector<Vec3>* positions, int d, int h, int a, int b ) {
    distHA.setPoints( positions, h, a );
    angleDHA.setPoints( positions, d, h, a );
    angleHAB.setPoints( positions, h, a, b );
  }
  virtual double energy();
  virtual void addToGradient_P(
    std::vector<Vec3> &dV_dr_P,
    std::vector<SecondDerivative> &secondDerivative_P );

  double calcAngleRadHAB() const { return acos(angleHAB.calcCosTheta()); }
  void setMinAngleRadHAB( double theta ) { angleHAB.setMinAngleRad( theta ); }
  void setkAngleHAB( double k ) { angleHAB.setk( k ); }
  int getb() const { return angleHAB.getp3(); }
  double getConstraintMinAngleRadHAB() const { return angleHAB.getMinAngleRad(); }
protected:
  MinAngleConstraint angleHAB;
};

#include "ConstraintContainer.h"
#include <limits>
class HBContainer : public ConstraintContainer<HBConstraint> {
public:
  HBContainer( int natoms ) : ConstraintContainer<HBConstraint>( natoms ) {}
  double worstDistanceViolation() const {
    double worst = 0;
    for ( const_iterator it = constraints.begin(); it != constraints.end(); it++ ) {
      double violation = it->calcDistHA() - it->getConstraintMaxDistHA();
      if ( violation > worst ) worst = violation;
      //if ( violation > 0.1 ) cout << "HB " << it->isBreakable() << " " << it->geth() << " "<< it->geta() << " " << violation << " " << it->getConstraintMaxDistHA() << endl;
    }
    return worst;
  }
  double worstAngleViolationDeg() const {
    double worst = 0;
    for ( const_iterator it = constraints.begin(); it != constraints.end(); it++ ) {
      double violation;
      if ( it->getConstraintMinAngleRadDHA() > std::numeric_limits<double>::epsilon() ) {
        violation = it->getConstraintMinAngleRadDHA() - it->calcAngleRadDHA();
        if ( violation > worst ) worst = violation;
      }
    }
    return worst*180.0/M_PI;
  }
  void write( std::string filename );
};

class BBHBContainer : public ConstraintContainer<BBHBConstraint> {
public:
  BBHBContainer( int natoms ) : ConstraintContainer<BBHBConstraint>( natoms ) {}
  double worstDistanceViolation() const {
    double worst = 0;
    for ( const_iterator it = constraints.begin(); it != constraints.end(); it++ ) {
      double violation = it->calcDistHA() - it->getConstraintMaxDistHA();
      if ( violation > worst ) worst = violation;
      //if ( violation > 0.1 ) cout << "BBHB " << it->isBreakable() << " " << it->geth() << " "<< it->geta() << " " << violation << " " << it->getConstraintMaxDistHA() << endl;
    }
    return worst;
  }
  double worstAngleViolationDeg() const {
    double worst = 0;
    for ( const_iterator it = constraints.begin(); it != constraints.end(); it++ ) {
      double violation;
      if ( it->getConstraintMinAngleRadDHA() > std::numeric_limits<double>::epsilon() ) {
        violation = it->getConstraintMinAngleRadDHA() - it->calcAngleRadDHA();
        if ( violation > worst ) worst = violation;
      }
      if ( it->getConstraintMinAngleRadHAB() > std::numeric_limits<double>::epsilon() ) {
        violation = it->getConstraintMinAngleRadHAB() - it->calcAngleRadHAB();
        if ( violation > worst ) worst = violation;
      }
    }
    return worst*180.0/M_PI;
  }
  void write( std::string filename );
};

#endif /* HBCONSTRAINT_H_ */
