// CorrPerturber.cpp
// Perturbs atoms along gradient of map correlation with respect to atom 
// position
// Craig Jolley,June 2008
////////////////////////////////////////////////////////////////////////////////

#include "RandomCenterPerturber.h"
#include "RandomVector.h"
#include "RigidUnitSystem.h"
#include "FitPerturber.h"
#include "CorrPerturber.h"
#include "Vec3.h"
#include <cmath>
#include "mt19937ar.h"
#include "GenericMap.h"

using namespace std;

CorrPerturber::CorrPerturber(RigidUnitSystem *rigidUnitSystem_,
			   GenericMap *map_, EMStructure *em_, double size_) {
  rigidUnitSystem = rigidUnitSystem_;
  map = map_;
  size = size_;
  fitPert = new FitPerturber(rigidUnitSystem);
  em = em_;
  symFitter = new Fit();
}

CorrPerturber::~CorrPerturber() {
  delete fitPert;
  delete symFitter;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Default perturbation routine; use this one when fitting without symmetry.
////////////////////////////////////////////////////////////////////////////////
void CorrPerturber::perturb() {
  vector <Vec3> newPositions = rigidUnitSystem->meanPositions();
  vector <Vec3> perturbations;
  perturbations.resize(newPositions.size());
  double maxGradNorm2 = 0;
  map->calcDiffMap(*em);
  //double total = 0;
  #pragma omp parallel
  {
    #pragma omp for
    for (size_t atomNum = 0; atomNum < newPositions.size(); atomNum++) {
      perturbations[atomNum] = map->cGradient(*em,atomNum);
      double norm2 = perturbations[atomNum].norm2(); 
      //total += sqrt(norm2);
      #pragma omp critical (MAP_PERTURB_MAX)
      {
        if (norm2 > maxGradNorm2) {
          maxGradNorm2 = norm2;
        }
      } 
    }
    #pragma omp single
    {
      if (maxGradNorm2 == 0) {
        cerr << "WARNING: No correlation gradient being felt by any atoms!" 
             << "  Skipping gradient perturbation.\n";
        perturbations.clear();
      }
    }
    #pragma omp for
    for (size_t atomNum = 0; atomNum < perturbations.size(); atomNum++) {
      newPositions[atomNum] += perturbations[atomNum] * (size / sqrt(maxGradNorm2));
    }
  } // end parallel region
  // atoms have been moved; use FitPerturber to fit rigid units
  fitPert->setMeanPointsTarget(&newPositions);
  fitPert->perturb();
  //cout << "Average map perturbation: " << total/sqrt(maxGradNorm2)/newPositions.size() << endl;
  return;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Perturbation routine for use with explicit symmetry
////////////////////////////////////////////////////////////////////////////////
void CorrPerturber::symPerturb() {
  vector <Vec3> newPositions = rigidUnitSystem->meanPositions();
  vector <Vec3> perturbations;
  perturbations.resize(newPositions.size());
  double maxGradNorm2 = 0;
  Vec3 netMotion(0,0,0);
  map->calcDiffMap(*em);
  //double total = 0;
  #pragma omp parallel
  {
    #pragma omp for
    for (size_t atomNum = 0; atomNum < newPositions.size(); atomNum++) {
      perturbations[atomNum] = map->cGradient(*em,atomNum);
      netMotion += perturbations[atomNum];
      double norm2 = perturbations[atomNum].norm2(); 
      //total += sqrt(norm2);
      #pragma omp critical (MAP_PERTURB_MAX)
      {
        if (norm2 > maxGradNorm2) {
          maxGradNorm2 = norm2;
        }
      } 
    }
    #pragma omp single
    {
      if (maxGradNorm2 == 0) {
        cerr << "WARNING: No correlation gradient being felt by any atoms!" 
             << "  Skipping gradient perturbation.\n";
        perturbations.clear();
      }
    }
    #pragma omp for
    for (size_t atomNum = 0; atomNum < perturbations.size(); atomNum++) {
      newPositions[atomNum] += perturbations[atomNum] * (size / sqrt(maxGradNorm2));
    }
  } // end parallel region
  // correct for possible rotation of symmetry axes
  symFitter->setTargetAbsolutePoints(newPositions);
  symFitter->setSourceAbsolutePoints(rigidUnitSystem->meanPositions());
  symFitter->simpleFit();
  //cout << "Rotor: " << symFitter->getFitRotor() << endl;
  Vec3 rotor = symFitter->getFitRotor();
  symMat->moveGlobalRotation(rotor);
  // correct for possible motion of symmetry origin
  netMotion *= (size / sqrt(maxGradNorm2));
  netMotion /= perturbations.size();
  symMat->moveOrigin(netMotion);
  //cout << "Average map perturbation: " << total/sqrt(maxGradNorm2)/newPositions.size() << endl;
  // atoms have been moved; use FitPerturber to fit rigid units
  fitPert->setMeanPointsTarget(&newPositions);
  fitPert->perturb();
  return;
}
