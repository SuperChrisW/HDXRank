/*
 * MinDistConstraintContainer.cpp
 *
 *  Created on: Mar 9, 2010
 *      Author: dwfarrel
 */

#include "MinDistConstraintContainer.h"
#include <iostream>
#include <fstream>

using namespace std;

double MinDistConstraintContainer::worstDistanceViolation() const {
  double worst = 0;
  for ( const_iterator it = constraints.begin(); it != constraints.end(); it++ ) {
    double violation = it->violation();
    if ( violation > worst ) worst = violation;
  }
  return worst;
}

void MinDistConstraintContainer::writeViolations() const {
  for ( const_iterator it = constraints.begin(); it != constraints.end(); it++ ) {
    double violation = it->violation();
    if ( violation > 0.02 ) //std::numeric_limits<double>::epsilon() )
      std::cout << it->getp1() << " " << it->getp2() << " " << violation << " " << it->getCutoff() << " " << it->calcDist() << std::endl;
  }
}

void MinDistConstraintContainer::write( string filename ) const {
  ofstream outfile( filename.c_str(), ios::out );
  for ( const_iterator it = constraints.begin(); it != constraints.end(); it++ ) {
    outfile << it->getp1() << " " << it->getp2() << " " << it->getCutoff() << '\n';
  }
  outfile.close();
}
void MinDistConstraintContainer::append( string filename ) const {
  ofstream outfile( filename.c_str(), ios::app );
  for ( const_iterator it = constraints.begin(); it != constraints.end(); it++ ) {
    outfile << it->getp1() << " " << it->getp2() << " " << it->getCutoff() << '\n';
  }
  outfile.close();
}
