#ifndef SYMMETRY_MATRICES_H_
#define SYMMETRY_MATRICES_H_

#include <vector>
#include <string>
#include "Vec3.h"
#include "Rotator.h"

typedef std::vector< std::vector<double> > Matrix;

class SymmetryMatrices {
public:
  SymmetryMatrices(std::string fileName);
  ~SymmetryMatrices();
  const Matrix getMatrix(int n) const;
  const Matrix getInverse(int n) const;
  const Matrix getProduct(int n, int m) const;
  const Vec3 getOffsetPair(int n, int m) const;
  size_t size() const;
  void moveOrigin(Vec3 &v) { origin += v; }
  const Vec3 getOrigin() const { return origin; }
  void moveGlobalRotation(Vec3 &v);
  const Matrix getGlobalRotation() const { return globalRotation; }
private:
  std::vector <Matrix> *matrices;
  std::vector <Matrix> *inverses;
  std::vector < std::vector <Matrix> > *products;
  // contains all pairwise products of matrices and inverses
  std::vector < std::vector <Vec3> > *offsetPairs;
  // contains all pairwise differences of offset vectors of matrices
  std::vector <Vec3> averagePos;
  // contains symmetry-averaged positions for all atoms
  Vec3 origin; 
  // origin of the coordinate system. (0,0,0) by default; can vary during
  // cryo-EM fitting
  Vec3 globalRotationRotor; // rotor containing global rotation
  Matrix globalRotation;
  Rotator rot1, rot2;
  // identity matrix by default; can vary during cryo-EM fitting
  Vec3 subtractOffsets(const Matrix &m1, const Matrix &m2);
  double determinant(const Matrix &m) const;
  void showMatrix(const Matrix &m) const;
  Matrix invert(const Matrix &m) const;
  Matrix mult(const Matrix &m1, const Matrix &m2) const;
  Matrix globalRotationTransform(const Matrix &m) const;
};

#endif /* SYMMETRY_MATRICES_H_ */
