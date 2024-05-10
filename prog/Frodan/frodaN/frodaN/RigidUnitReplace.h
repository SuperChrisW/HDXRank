/*
 * RigidUnitReplace.h
 *
 *  Created on: Sep 12, 2008
 *      Author: dwfarrel
 */

#ifndef RIGIDUNITREPLACE_H_
#define RIGIDUNITREPLACE_H_

#include "Fit.h"
class RigidUnitSystem;
class TargetEnergy;
#include <vector>
#include <queue>
class NeighborTable;
#include <map>

class RigidUnitReplace {
public:
  RigidUnitReplace(
      RigidUnitSystem *rigidUnitSystem,
      const TargetEnergy *targetEnergy,
      const std::map<int,int> *atommap,
      const NeighborTable* covTableA,
      const NeighborTable* covTableB );
  virtual ~RigidUnitReplace();
  void replace( int& count );
  void initializePQ( double rmsdupper, double rmsdlower );
  bool isReplaceableRU( int ru ) const { return isReplaceable[ru]; }
  void replaceRU( int ru );
private:
  bool checkRigidUnit( int ru,
      const std::map<int,int> *atommap,
      const NeighborTable* covTableA,
      const NeighborTable* covTableB );
  RigidUnitSystem *rigidUnitSystem;
  const TargetEnergy *targetEnergy;
  Fit fit;
  class PQNode;
  std::priority_queue<PQNode> pq;
  std::vector<char> isReplaceable;
};

class RigidUnitReplace::PQNode {
public:
  double rmsd;
  int ru;
  bool operator<( const RigidUnitReplace::PQNode& other ) const {
    return rmsd < other.rmsd;
  }
};



#endif /* RIGIDUNITREPLACE_H_ */
