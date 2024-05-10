/*
 * Repulsion.h
 *
 *  Created on: Jul 27, 2009
 *      Author: dwfarrel
 */

#ifndef REPULSION_H_
#define REPULSION_H_

#include <set>
#include <vector>
class RigidUnitSystem;
class NeighborTable;
class PairTypeCutoffs;
#include "TestForExclusion.h"
#include "ProteinInfo.h"
#include "Vec3.h"

class Repulsion {
public:
  Repulsion(
    const RigidUnitSystem* sys,
    const NeighborTable* nt,
    PairTypeCutoffs* pairTypeCutoffs );
  virtual ~Repulsion();

  void getCutoff( int p1, int p2, bool& isPairExcluded, double& cutoff ) const;
  double getMaxCutoff() const { return maxCutoff; }
  void exclude( int i, int j );

private:
  TestForExclusion testForExclusion;
  PairTypeCutoffs* pairTypeCutoffs;
  double maxCutoff;
  std::vector< std::set<int> > lookupExcludePair_I_LT_J;

};

#endif /* REPULSION_H_ */
