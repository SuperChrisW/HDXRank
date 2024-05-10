// MapPerturber.h
// Perturbs atoms along calculated gradient of electron density map
// Craig Jolley, February 2008
////////////////////////////////////////////////////////////////////////////////

#ifndef MAP_PERTURBER_H_
#define MAP_PERTURBER_H_

#include "RigidUnitSystem.h"
#include "FitPerturber.h"
#include "Vec3.h"
#include "GenericMap.h"
#include "SymmetryMatrices.h"
#include "Fit.h"
#include "Rotator.h"

class RigidUnitSystem;

class MapPerturber {
public:
  MapPerturber(RigidUnitSystem *rigidUnitSystem_, 
	      GenericMap *map_, double size_);
  virtual ~MapPerturber();
  void perturb();
  void associateMatrices(SymmetryMatrices *symMat_) {symMat = symMat_;}
  void symPerturb();
private:
  SymmetryMatrices *symMat;
  RigidUnitSystem *rigidUnitSystem;
  GenericMap *map;
  double size;
  FitPerturber *fitPert;
  Fit *symFitter;
};

#endif /* MAP_PERTURBER_H_ */
