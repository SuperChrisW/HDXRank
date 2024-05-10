#ifndef GRID_H_
#define GRID_H_

#include <vector>
#include "RigidUnitSystem.h"
#include "Vec3.h"
#include "PDB.h"

class Grid {
public:
  Grid(RigidUnitSystem &rigidUnitSystem_, PDB &pdb, double gridLength_);
  void updateGrid();
  std::vector <int> nearbyAtoms(Vec3 v);
private:
  RigidUnitSystem *rigidUnitSystem;
  double gridLength;
  Vec3 min, max;
  int sizeX, sizeY, sizeZ;
  void findMinMax();
  void setupGrid();
  int VecToArrayIndex(Vec3 v);
  std::vector <int> contributingAtoms;
  std::vector <std::vector <int> > myGrid;
  std::vector <std::vector <int> > adjacentPoints;
};

#endif
