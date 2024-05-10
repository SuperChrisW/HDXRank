// MapPerturber.cpp
// Perturbs atoms along calculated gradient of electron density map
// Craig Jolley, February 2008
////////////////////////////////////////////////////////////////////////////////

#include "RandomCenterPerturber.h"
#include "RandomVector.h"
#include "RigidUnitSystem.h"
#include "FitPerturber.h"
#include "MapPerturber.h"
#include "Vec3.h"
#include <cmath>
#include "mt19937ar.h"
#include "GenericMap.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Constructor for MapPerturber class
// Parameters:
//   RigidUnitSystem *rigidUnitSystem_ -- RigidUnitSystem to be perturbed
//   GenericMap *map_ -- EM map used to guide perturbations
//   double size_ -- perturbation size
////////////////////////////////////////////////////////////////////////////////
MapPerturber::MapPerturber(RigidUnitSystem *rigidUnitSystem_,
                         GenericMap *map_,double size_) {
  rigidUnitSystem = rigidUnitSystem_;
  map = map_;
  size = size_;
  fitPert = new FitPerturber(rigidUnitSystem);
  symFitter = new Fit();
}


////////////////////////////////////////////////////////////////////////////////
// Description:
//   Destructor for MapPerturber class.  Just frees up some memory allocated
//   by the constructor.
////////////////////////////////////////////////////////////////////////////////
MapPerturber::~MapPerturber() {
  delete fitPert;
  delete symFitter;
}


////////////////////////////////////////////////////////////////////////////////
// Description:
//   Default perturbation routine; use this one when fitting without symmetry.
////////////////////////////////////////////////////////////////////////////////
void MapPerturber::perturb() {
  vector <Vec3> newPositions = rigidUnitSystem->meanPositions();
  vector <Vec3> perturbations;
  perturbations.resize(newPositions.size());
  double maxGradNorm2 = 0;
  //double total = 0;
  #pragma omp parallel
  {
    #pragma omp for
    for (size_t atomNum = 0; atomNum < newPositions.size(); atomNum++) {
      perturbations[atomNum] = map->mapGradient(newPositions[atomNum]);
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
        cerr << "WARNING: No map gradient being felt by any atoms!  Skipping gradient perturbation.\n";
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
//   Perturbation routine for use with symmetry.  
////////////////////////////////////////////////////////////////////////////////
void MapPerturber::symPerturb() {
  vector <Vec3> newPositions = rigidUnitSystem->meanPositions();
  vector <Vec3> perturbations;
  perturbations.resize(newPositions.size());
  double maxGradNorm2 = 0;
  size_t maxGradAt;
  Vec3 netTranslation(0,0,0);
  //double total = 0;
  #pragma omp parallel
  {
    #pragma omp for
    for (size_t atomNum = 0; atomNum < newPositions.size(); atomNum++) {
      perturbations[atomNum] = map->mapGradient(newPositions[atomNum]);
      netTranslation += perturbations[atomNum];
      double norm2 = perturbations[atomNum].norm2(); 
      //total += sqrt(norm2);
      #pragma omp critical (MAP_PERTURB_MAX)
      {
        if (norm2 > maxGradNorm2) {
          maxGradNorm2 = norm2;
          maxGradAt = atomNum;
        }
      } 
    }
    #pragma omp single
    {
      if (maxGradNorm2 == 0) {
        cerr << "WARNING: No map gradient being felt by any atoms!  Skipping gradient perturbation.\n";
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
  // correct for possible translation of symmetry origin
  netTranslation *= (size / sqrt(maxGradNorm2));
  netTranslation /= perturbations.size();
  //cout << "New origin: " << symMat->getOrigin() << endl;
  symMat->moveOrigin(netTranslation);
  // atoms have been moved; use FitPerturber to fit rigid units
  fitPert->setMeanPointsTarget(&newPositions);
  fitPert->perturb();
  //cout << "Average map perturbation: " << total/sqrt(maxGradNorm2)/newPositions.size() << endl;
  //cout << "New origin: " << symMat->getOrigin() << endl;
  return;
}
