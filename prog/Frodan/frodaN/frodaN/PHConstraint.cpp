/*
 * PHConstraint.cpp
 *
 *  Created on: Sep 23, 2009
 *      Author: dwfarrel
 */

#include "PHConstraint.h"
#include <fstream>
#include <string>

using namespace std;

void PHContainer::write( string filename ) {
  ofstream outfile( filename.c_str(), ios::out );
  for ( const_iterator it = constraints.begin(); it != constraints.end(); it++ ) {
    outfile << it->getp1() << " " << it->getp2() << " " << it->getCutoff() << '\n';
  }
  outfile.close();
}

