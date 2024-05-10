/*
 * MinDistConstraintContainer.h
 *
 *  Created on: Mar 9, 2010
 *      Author: dwfarrel
 */

#ifndef MINDISTCONSTRAINTCONTAINER_H_
#define MINDISTCONSTRAINTCONTAINER_H_

#include "ConstraintContainer.h"
#include "DistConstraint.h"

class MinDistConstraintContainer : public ConstraintContainer<MinDistConstraint> {
public:
  MinDistConstraintContainer( int natoms ) : ConstraintContainer<MinDistConstraint>( natoms ) {}
  virtual ~MinDistConstraintContainer() {}
  double worstDistanceViolation() const;
  void writeViolations() const;
  void write( std::string filename ) const;
  void append( std::string filename ) const;

};

#endif /* MINDISTCONSTRAINTCONTAINER_H_ */
