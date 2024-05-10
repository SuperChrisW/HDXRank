/*
 * CollectiveRandomPerturber.h
 *
 *  Created on: Apr 27, 2010
 *      Author: dwfarrel
 */

#ifndef COLLECTIVERANDOMPERTURBER_H_
#define COLLECTIVERANDOMPERTURBER_H_

#include <vector>
#include <set>
class ProteinInfo;
#include "NeighborTable.h"
class RigidUnitSystem;
class ConstraintEnforcingPotential;

class CollectiveRandomPerturber {
public:
  CollectiveRandomPerturber(
    const ProteinInfo* prot_,
    const NeighborTable* covNT_,
    RigidUnitSystem* sys_,
    const ConstraintEnforcingPotential* cep_ );
  virtual ~CollectiveRandomPerturber();
  void perturb();

private:
  const ProteinInfo* prot;
  const NeighborTable& covNT;
  RigidUnitSystem* sys;
  const ConstraintEnforcingPotential* cep;
  NeighborTable resNT;
  std::vector< std::vector<int> > RUlistFromResi;
  std::set<int> resi_bucket;

  class PQNode;

  void perturbRigidUnitGroup_Translational( const std::set<int>& rigidUnitGroup );
  void perturbRigidUnitGroup_Rotational( const std::set<int>& rigidUnitGroup );
  void extendGroup( int resi, int nExtend, std::vector<int>& resi_group );

};

class CollectiveRandomPerturber::PQNode {
public:
  PQNode() {}
  PQNode( double p, int r ) : priority(p), resi(r) {}
  double priority;
  int resi;
  bool operator<( const CollectiveRandomPerturber::PQNode& other ) const {
    return priority < other.priority;
  }
};



#endif /* COLLECTIVERANDOMPERTURBER_H_ */
