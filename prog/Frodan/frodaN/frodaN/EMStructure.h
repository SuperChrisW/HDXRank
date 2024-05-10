// Author name:  Craig Jolley
// Created:      02 Jun 2008

#ifndef EM_STRUCTURE_H_
#define EM_STRUCTURE_H_

#include <vector>
#include "Vec3.h"
#include "PDB.h"
#include "RigidUnitSystem.h"
#include "Grid.h"

//////////////////
// This code was originally written for FRODA (beginning in June 2006) and is 
// being adapted for cryo-EM fitting using newsim.
// Craig Jolley, June 2008
//////////////////

class EMStructure {
private: 
  struct Atom {
    int id; // ID in rigidUnitSystem
    short int Z;
  };
private:
  double distanceSquared(int id, Vec3 pos);
  Grid *grid;
  RigidUnitSystem *rigidUnitSystem;
  double sigma; // resolution / 2
  std::vector<double> lookupTable;
  double binSize;
  std::vector<Atom> atoms; 
public:
  double cutoff;  
  EMStructure(PDB &pdb, RigidUnitSystem &RUSInput, double resolution, double contained);
  ~EMStructure();
  double scaledDensity(double r2, int atomID);
  double densityAt(Vec3 v);
  void updateGrid() { grid->updateGrid(); }  
  Vec3 atomPos(int id) { return rigidUnitSystem->meanPositions(id); }
};

#endif
