/*
 * HBConstraint.cpp
 *
 *  Created on: Jul 20, 2009
 *      Author: dwfarrel
 */

#include "HBConstraint.h"
#include <limits>
#include <fstream>
#include <string>

using namespace std;

//Abstract Class
HBConstraintAbstract::HBConstraintAbstract() :
  isBreakable_( true )
{
}

HBConstraintAbstract::~HBConstraintAbstract() {}

//HB Class
HBConstraint::HBConstraint() {}

HBConstraint::HBConstraint( const std::vector<Vec3>* positions, int d, int h, int a ) {
  setAtoms( positions, d, h, a );
}

HBConstraint::~HBConstraint() {}

double HBConstraint::energy() {
  return distHA.energy() +
         ( angleDHA.getMinAngleRad() > numeric_limits<double>::epsilon() ? angleDHA.energy() : 0 );
}

void HBConstraint::addToGradient_P(
  std::vector<Vec3> &dV_dr_P,
  std::vector<SecondDerivative> &secondDerivative_P ) {

  distHA.addToGradient_P( dV_dr_P, secondDerivative_P );
  if ( angleDHA.getMinAngleRad() > numeric_limits<double>::epsilon() ) angleDHA.addToGradient_P( dV_dr_P, secondDerivative_P );
}

//BBHB Class
BBHBConstraint::BBHBConstraint() {}

BBHBConstraint::~BBHBConstraint() {}

BBHBConstraint::BBHBConstraint( const std::vector<Vec3>* positions, int d, int h, int a, int b ) {
  setAtoms( positions, d, h, a, b );
}

double BBHBConstraint::energy() {
  return distHA.energy() +
         ( angleDHA.getMinAngleRad() > numeric_limits<double>::epsilon() ? angleDHA.energy() : 0 ) +
         ( angleHAB.getMinAngleRad() > numeric_limits<double>::epsilon() ? angleHAB.energy() : 0 );
}

void BBHBConstraint::addToGradient_P(
  std::vector<Vec3> &dV_dr_P,
  std::vector<SecondDerivative> &secondDerivative_P ) {

  distHA.addToGradient_P( dV_dr_P, secondDerivative_P );
  if ( angleDHA.getMinAngleRad() > numeric_limits<double>::epsilon() ) angleDHA.addToGradient_P( dV_dr_P, secondDerivative_P );
  if ( angleHAB.getMinAngleRad() > numeric_limits<double>::epsilon() ) angleHAB.addToGradient_P( dV_dr_P, secondDerivative_P );
}

void HBContainer::write( string filename ) {
  ofstream outfile( filename.c_str(), ios::out );
  for ( const_iterator it = constraints.begin(); it != constraints.end(); it++ ) {
    outfile << it->geth() << " " << it->geta() << " " << it->getConstraintMaxDistHA() << " " << it->getConstraintMinAngleRadDHA() << '\n';
  }
  outfile.close();
}

void BBHBContainer::write( string filename ) {
  ofstream outfile( filename.c_str(), ios::out );
  for ( const_iterator it = constraints.begin(); it != constraints.end(); it++ ) {
    outfile << it->geth() << " " << it->geta() << " " << it->getConstraintMaxDistHA() << " "
            << it->getConstraintMinAngleRadDHA() << " " << it->getConstraintMinAngleRadHAB() << '\n';
  }
  outfile.close();
}
