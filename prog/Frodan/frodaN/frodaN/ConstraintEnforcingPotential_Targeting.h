/*
 * TargetingConstraintManager.h
 *
 *  Created on: Sep 17, 2009
 *      Author: dwfarrel
 */

#ifndef CONSTRAINTENFORCINGPOTENTIAL_TARGETING_H_
#define CONSTRAINTENFORCINGPOTENTIAL_TARGETING_H_

#include "ConstraintEnforcingPotential.h"
#include "TargetMap.h"
#include "NeighborTable.h"
#include "TargetEnergy.h"
#include <vector>
#include <map>
#include <set>
#include <string>
#include "Vec3.h"
class Settings;
class RigidUnitSystem;
class ProteinInfo;

class ConstraintEnforcingPotential_Targeting : public ConstraintEnforcingPotential {
public:
  ConstraintEnforcingPotential_Targeting(
    RigidUnitSystem* sys,
    const ProteinInfo& prot,
    const NeighborTable& nt,
    const Settings& settings );

  virtual ~ConstraintEnforcingPotential_Targeting();

  virtual std::string generateColumnHeaderString();
  virtual std::string generateSummaryString();

  TargetEnergy *targetEnergy;
  TargetMap* targmap;
  NeighborTable* targneigh;

private:
  void setupConstraintsBBHB( BBHBContainer* bbhb, BBHBContainer* targbbhb );
  void setupConstraintsHB( HBContainer* hb, HBContainer* targhb );
  void setupConstraintsPH( PHContainer* ph, PHContainer* targph );
  void setupConstraintsSC( SideChainTorsionContainer* tor, SideChainTorsionContainer* targtor );
  void addMinDistOverridesFromTarg( const ConstraintEnforcingPotential* targcep );
  void initialize(
    RigidUnitSystem* sys,
    const ProteinInfo& prot,
    const Settings& settings,
    const std::vector<Vec3>& targcoords );
  bool fixedCommon;
  bool fixedNoncommon;
  bool removeCommon;
  bool removeNoncommon;
  std::set< int > breakable_atomlist;

};

#endif /* CONSTRAINTENFORCINGPOTENTIAL_TARGETING_H_ */
