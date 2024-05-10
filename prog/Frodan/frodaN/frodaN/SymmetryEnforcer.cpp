// SymmetryEnforcer.cpp
// Energy term to enforce symmetry
// Craig Jolley, January 2008
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include "SymmetryEnforcer.h"
#include "RigidUnitSystem.h"
#include <iomanip>
#include <limits>

using namespace std;

SymmetryEnforcer::SymmetryEnforcer(const RigidUnitSystem *rigidUnitSystem_,
                 const SymmetryMatrices *symmetryMatrices_) {
  rigidUnitSystem = rigidUnitSystem_;
  symmetryMatrices = symmetryMatrices_;
  averagePositions.resize(rigidUnitSystem->nPoints());
  updateAvgPositions();
}  
          
SymmetryEnforcer::~SymmetryEnforcer() {}

Vec3 SymmetryEnforcer::transform(Vec3 v, Matrix m) {
  // this procedure doesn't use offset vectors; those need
  // to be dealt with explicitly
  Vec3 newV;
  v -= symmetryMatrices->getOrigin();
  // correct for possible global rotations during cryo-EM fitting
  newV.x = m[0][0]*v.x + m[0][1]*v.y + m[0][2]*v.z;
  newV.y = m[1][0]*v.x + m[1][1]*v.y + m[1][2]*v.z;
  newV.z = m[2][0]*v.x + m[2][1]*v.y + m[2][2]*v.z;
  newV += symmetryMatrices->getOrigin();
  return newV;
}

Vec3 SymmetryEnforcer::getAvgPosition(int p) {
  Vec3 result;
  result.x = result.y = result.z = 0;
  size_t nMatrices = symmetryMatrices->size();
  int monomerAtoms = rigidUnitSystem->nPoints() / nMatrices;
  int home = (int) p / monomerAtoms;
  int index = p % monomerAtoms;
  for (size_t j = 0; j < nMatrices; j++) {
    const Vec3 neighborPos = rigidUnitSystem->meanPositions(j*monomerAtoms + index);
    result += transform(neighborPos,symmetryMatrices->getProduct(home,j));
    result += symmetryMatrices->getOffsetPair(home,j);
  }
  result /= nMatrices;
  return result;
}

double SymmetryEnforcer::energy() {
  double totalEnergy = 0.0;
  for (size_t p = 0; p < rigidUnitSystem->nPoints(); p++) {
    Vec3 meanPos = rigidUnitSystem->meanPositions(p);
    totalEnergy += meanPos.dist2(averagePositions[p]);
  }
  totalEnergy /= 2; // because energy is (1/2)kx^2
  //cout << "Symmetry energy: " << totalEnergy << endl;
  return totalEnergy;
}

void SymmetryEnforcer::addToGradient(std::vector<Vec3> &dV_dr_rigidUnitPoint,
				     std::vector<SecondDerivative> &secondDerivative_rigidUnitPoint) {
  // For details of the math here, see Craig's notebook entries for 01/10/08 - 01/21/08
  size_t nMatrices = symmetryMatrices->size();
  int monomerAtoms = rigidUnitSystem->nPoints() / nMatrices;
  // diagonal elements of second-derivative matrix are all equal
  double diagElement = 1.0 - 1.0 / nMatrices;
  #pragma omp parallel
  {
    #pragma omp for 
    for (size_t p = 0; p < rigidUnitSystem->nPoints(); p++) {
      // gradient vector is based on RUPs, but all RUPs associated
      // with a given point will have the same contribution to the gradient
      int home = (int) p / monomerAtoms;
      int index = p % monomerAtoms;
      Vec3 delta(0,0,0); // this will be added to first derivative
      // only gradient terms for an RUP and its symmetry relatives are nonzero
      for (size_t j = 0; j < nMatrices; j++) {
        int replicaP = index + j*nMatrices;
        Vec3 v1 = averagePositions[replicaP] - rigidUnitSystem->meanPositions(replicaP);
        // correction for global rotations of symmetry axes done within symmetryMatrices
        Matrix m = symmetryMatrices->getProduct(j,home);
        // get x derivative
        Vec3 v2(m[0][0],m[1][0],m[2][0]);
        delta.x += v1.dot(v2);
        // now y derivative
        v2.x = m[0][1]; v2.y = m[1][1]; v2.z = m[2][1];
        delta.y += v1.dot(v2);
        // and z derivative
        v2.x = m[0][2]; v2.y = m[1][2]; v2.z = m[2][2];
        delta.z += v1.dot(v2);    
      }
      delta /= nMatrices;
      delta -= averagePositions[p] - rigidUnitSystem->meanPositions(p);
      vector <int> rups = rigidUnitSystem->getRUPlistFromP(p);
      delta /= rups.size();
      double secondD = diagElement / (rups.size()*rups.size());
      // TODO: check lab notebook; why did I have this squared again?
      for (size_t i = 0; i < rups.size(); i++) {
        dV_dr_rigidUnitPoint[rups[i]] += delta;
        secondDerivative_rigidUnitPoint[rups[i]].d2V_dx2 = secondD;
        secondDerivative_rigidUnitPoint[rups[i]].d2V_dy2 = secondD;
        secondDerivative_rigidUnitPoint[rups[i]].d2V_dz2 = secondD;
      }
    }
  } // end parallel region
  return;
}

double SymmetryEnforcer::mismatch() {
  double maxMismatch2 = 0;
  double dist2;
  size_t savedP = 0;
  for (size_t p = 0; p < rigidUnitSystem->nPoints(); p++) {
    Vec3 meanPos = rigidUnitSystem->meanPositions(p);
    dist2 = meanPos.dist2(averagePositions[p]);
    if (dist2 > maxMismatch2) {
      maxMismatch2 = dist2;
      savedP = p;
    }
  }
  if (verbose && maxMismatch2 > numeric_limits<double>::epsilon()) {
    cout << "Max distance to symmetric average: atom " <<
      savedP << " dist " << sqrt(maxMismatch2) << endl;
  }
  return sqrt(maxMismatch2);
}

void SymmetryEnforcer::updateAvgPositions() {
  size_t nMatrices = symmetryMatrices->size();
  size_t monomerAtoms = rigidUnitSystem->nPoints() / nMatrices;
  #pragma omp parallel for
  for (size_t p = 0; p < monomerAtoms; p++) {
    averagePositions[p] = getAvgPosition(p);
    for (size_t m = 1; m < nMatrices; m++) {
      averagePositions[m*monomerAtoms + p] = transform(averagePositions[p],
                                                    symmetryMatrices->getMatrix(m));
    }
  }
}

