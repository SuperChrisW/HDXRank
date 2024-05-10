/*
 * SideChainTorsion.h
 *
 *  Created on: Nov 14, 2009
 *      Author: dwfarrel
 */

#ifndef SIDECHAINTORSION_H_
#define SIDECHAINTORSION_H_

#include "MinDistConstraintContainer.h"
#include <vector>
class NeighborTable;
class ProteinInfo;
#include "Vec3.h"
#include <set>
#include <string>
class RigidUnitSystem;

typedef MinDistConstraintContainer SideChainTorsionContainer;

class SideChainTorsionInitializer {
public:
  SideChainTorsionInitializer(
    const ProteinInfo& prot,
    const NeighborTable& nt,
    const std::vector<Vec3>& coords,
    SideChainTorsionContainer& tor );
  virtual ~SideChainTorsionInitializer();
  void setupConstraints();
  void excludeRigidPairs( const RigidUnitSystem* sys );
private:
  const ProteinInfo& prot;
  const NeighborTable& nt;
  const std::vector<Vec3>& coords;
  SideChainTorsionContainer& tor;
  double k;
  std::set<std::string> resList;

  double calc14Constraint(
    int a1, int a2, int a3, int a4,
    double costheta, double sintheta );
};

#endif /* SIDECHAINTORSION_H_ */
