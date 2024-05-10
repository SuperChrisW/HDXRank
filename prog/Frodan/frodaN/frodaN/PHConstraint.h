/*
 * PHConstraint.h
 *
 *  Created on: Aug 24, 2009
 *      Author: dwfarrel
 */

#ifndef PHCONSTRAINT_H_
#define PHCONSTRAINT_H_
#include <iostream>

using namespace std;

#include "DistConstraint.h"

class PHConstraint : public MaxDistConstraint {
public:
  PHConstraint() : isBreakable_( true ) {}
  virtual ~PHConstraint() {}
  void makeUnbreakable() { isBreakable_ = false; }
  bool isBreakable() const { return isBreakable_; }
private:
  bool isBreakable_;
};

#include "ConstraintContainer.h"
class PHContainer : public ConstraintContainer<PHConstraint> {
public:
  PHContainer( int natoms ) : ConstraintContainer<PHConstraint>( natoms ) {}
  double worstDistanceViolation() const {
    double worst = 0;
    for ( const_iterator it = constraints.begin(); it != constraints.end(); it++ ) {
      double violation = it->violation();
      if ( violation > worst ) worst = violation;
      //if ( violation > 0.1 ) cout << "PH " << it->isBreakable() << " " << it->getp1() << " "<< it->getp2() << " " << violation << " " << it->getCutoff() << endl;
    }
    return worst;
  }

  void write( std::string filename );
};

#endif /* PHCONSTRAINT_H_ */
