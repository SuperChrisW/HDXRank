#ifndef RIGIDUNITDATA_H_
#define RIGIDUNITDATA_H_

#include "Vec3.h"
#include <vector>
#include "AbstractRigidUnitData.h"

class RigidUnitData : public AbstractRigidUnitData
{
public:
  RigidUnitData() {}
  virtual ~RigidUnitData() {}

  void setSizes( int nP, int nRU, int nRUP ) {
    rotors.resize( nRU );
    centers.resize( nRU );
    basePositions.resize( nRUP );
    radius.resize( nRU );
    hasZeroRadius.resize( nRU );
    absolutePositions.resize( nRUP );
    meanPositions.resize( nP );
  }

  std::vector<Vec3> meanPositions;

};

#endif /*RIGIDUNITDATA_H_*/
