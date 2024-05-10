/*
 * RamaOverlapEnergy.h
 *
 *  Created on: Nov 13, 2009
 *      Author: dwfarrel
 */

#ifndef RAMA_H_
#define RAMA_H_

#include <vector>
class ProteinInfo;
class NeighborTable;
#include "Vec3.h"
#include "ConstraintContainer.h"
#include "DistConstraint.h"
class RigidUnitSystem;
#include <set>
#include <string>

class RamaContainer : public ConstraintContainer<MinDistConstraint> {
public:

  RamaContainer(
    const RigidUnitSystem *sys,
    const ProteinInfo& proteinInfo,
    const NeighborTable& nt );
  virtual ~RamaContainer();

  double worstDistanceViolation() const {
    double worst = 0;
    for ( const_iterator it = constraints.begin(); it != constraints.end(); it++ ) {
      double violation = it->getCutoff() - it->calcDist();
      if ( violation > worst ) worst = violation;
      //if ( violation > 0.1 ) cout << "SC " << it->getp1() << " "<< it->getp2() << " " << violation << " " << it->getCutoff() << endl;
    }
    return worst;
  }
  void writeViolations() const {
    std::cout << "Rama Violations: " << std::endl;
    for ( const_iterator it = constraints.begin(); it != constraints.end(); it++ ) {
      double violation = it->getCutoff() - it->calcDist();
      if ( violation > 0.02 ) //std::numeric_limits<double>::epsilon() )
        std::cout << it->getp1() << " " << it->getp2() << " " << violation << " " << it->getCutoff() << " " << it->calcDist() << std::endl;
    }
  }


private:
  const RigidUnitSystem* sys;
  const ProteinInfo& proteinInfo;
  const NeighborTable& nt;
  double k;
  std::set<std::string> resList;
  void attemptAddRamaConstraint( int i, int j, double c );
};

#endif /* RAMA_H_ */
