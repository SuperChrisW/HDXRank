// CorrPerturber.h
// Perturbs atoms along gradient of correlation coefficient with respect
// to atom positions.
// Craig Jolley, June 2008
////////////////////////////////////////////////////////////////////////////////

#ifndef CORR_PERTURBER_H_
#define CORR_PERTURBER_H_

#include "RigidUnitSystem.h"
#include "FitPerturber.h"
#include "Vec3.h"
#include "GenericMap.h"
#include "EMStructure.h"
#include "SymmetryMatrices.h"

class RigidUnitSystem;

class CorrPerturber {
public:
  CorrPerturber(RigidUnitSystem *rigidUnitSystem_, 
		GenericMap *map_, EMStructure *em_, double size_);
  virtual ~CorrPerturber();
  void perturb();
  void associateMatrices(SymmetryMatrices *symMat_) {symMat = symMat_;}
  void symPerturb();
private:
  SymmetryMatrices *symMat;
  RigidUnitSystem *rigidUnitSystem;
  GenericMap *map;
  double size;
  FitPerturber *fitPert;
  EMStructure *em;
  Fit *symFitter;
};

#endif /* MAP_PERTURBER_H_ */
