// SymmetryMatrices.cpp
// Stores matrices, inverses, and their products
// Craig Jolley, January 2008
////////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>
#include "SymmetryMatrices.h"

using namespace std;

SymmetryMatrices::SymmetryMatrices(string fileName) {
  // Read in matrices from file
  // NOTE:  The "BIOMT" format isn't necessarily standard; this may require 
  // adaptation in some cases
  ifstream inFile;
  inFile.open(fileName.c_str());
  if (!inFile.is_open()) {
    cerr << "Could not open " << fileName << " to read BioMT matrices!\n";
    exit(EXIT_FAILURE);
  }
  matrices = new vector<Matrix>;
  while (inFile.good()) {
    string input;  
    getline(inFile,input);
    if (input.substr(0,19) == "REMARK 350   BIOMT1") {
      vector <double>  nullVec(4,0.0);
      Matrix tempMatrix(3,nullVec);
      for (int i = 0; i < 3; i++) {
        tempMatrix[0][i] = atof(input.substr(24+i*10,9).c_str());
      }
      tempMatrix[0][3] = atof(input.substr(61,7).c_str());
      getline(inFile,input); // get second line
      if (input.substr(0,19) != "REMARK 350   BIOMT2") {
        cerr << "ERROR: BIOMT1 should be followed by BIOMT2! Skipping matrix.\n";
      } else {
        for (int i = 0; i < 3; i++) {
          tempMatrix[1][i] = atof(input.substr(24+i*10,9).c_str());
        }
        tempMatrix[1][3] = atof(input.substr(61,7).c_str());
        getline(inFile,input); // get third line
        if (input.substr(0,19) != "REMARK 350   BIOMT3") {
          cerr << "ERROR: BIOMT2 should be followed by BIOMT3! Skipping matrix.\n";
        } else {
          for (int i = 0; i < 3; i++) {
            tempMatrix[2][i] = atof(input.substr(24+i*10,9).c_str());
          }
          tempMatrix[2][3] = atof(input.substr(60,8).c_str());
          // tempMatrix is now complete; save it
          matrices->push_back(tempMatrix);
        }
      }
    }
  }
  if (!inFile.eof() && inFile.fail()) {
    cerr << "Error in reading file " << fileName << "; BIOMT matrices may not be complete.\n";
  }
  inFile.close();
  cout << "Read " << matrices->size() << " matrices from " << fileName << endl;
  // if matrices[0] isn't the identity matrix, there will be trouble
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 4; j++) {
      if (i == j && matrices->at(0)[i][j] != 1.0) {
        cerr << "ERROR: First BIOMT entry should be the identity matrix!\n";
        exit(EXIT_FAILURE);
      } else if (i != j && matrices->at(0)[i][j] != 0.0) {
        cerr << "ERROR: First BIOMT entry should be the identity matrix!\n";
        exit(EXIT_FAILURE);
      }        
    }    
  }
  // Matrices have been loaded; now generate inverses and products
  inverses = new vector <Matrix>;
  products = new vector < vector <Matrix> >;
  offsetPairs = new vector < vector <Vec3> >;
  for (size_t i = 0; i < matrices->size(); i++) {
    inverses->push_back(invert(matrices->at(i)));
  }
  products->resize(matrices->size());
  offsetPairs->resize(matrices->size());
  for (size_t i = 0; i < matrices->size(); i++) {
    products->at(i).resize(matrices->size());
    offsetPairs->at(i).resize(matrices->size());
    for (size_t j = 0; j < matrices->size(); j++) {
      products->at(i).at(j) = mult(matrices->at(i),inverses->at(j));
      offsetPairs->at(i).at(j) = subtractOffsets(matrices->at(i),matrices->at(j));
    }
  }
  origin = Vec3(0,0,0); // this will remain zero unless doing cryo-EM fitting
  globalRotation.resize(3);
  for (int i = 0; i < 3; i++) {
    globalRotation[i].resize(3,0.0);
    globalRotation[i][i] = 1.0;
  } // globalRotation is now the identity matrix; may change during cryo-EM fitting
}

SymmetryMatrices::~SymmetryMatrices() {
  delete inverses;
  delete products;
  delete offsetPairs;
  delete matrices;
}


void SymmetryMatrices::moveGlobalRotation(Vec3 &v) {
  rot1.setRotor(globalRotationRotor);
  rot2.setRotor(v);
  rot1.addRotation(rot2);
  globalRotationRotor = rot1.rotor();
  globalRotation = rot1.getRotationMatrix();
}


Matrix SymmetryMatrices::mult(const Matrix &m1, const Matrix &m2) const {
  // For calculating gradients I only need the rotational parts of the
  // matrices, not the constant offset vectors
  Matrix result; // make this 3x3
  result.resize(3);
  for (int i = 0; i < 3; i++) {
    result[i].resize(3);
    for (int j = 0; j < 3; j++) {
      result[i][j] = 0;
      for (int k = 0; k < 3; k++) {
        result[i][j] += m1[i][k]*m2[k][j];
      }
    }
  }
  return result;
}

Vec3 SymmetryMatrices::subtractOffsets(const Matrix &m1, const Matrix &m2) {
  // if o(n) is the vector offset and m(n) is the rotation matrix, I want
  // o(m1) - m(m1)o(m2)
  Vec3 result;
  result.x = m1[0][3];
  result.y = m1[1][3];
  result.z = m1[2][3];
  for (int k = 0; k < 3; k++) {
    result.x -= m1[k][0]*m2[0][3];
    result.y -= m1[k][1]*m2[1][3];
    result.z -= m1[k][2]*m2[2][3];
  }
  return result;
}


double SymmetryMatrices::determinant(const Matrix &m) const {
  // Assumes it's being given a 3x3 matrix
  double result = 0;
  result += m[0][0]*m[1][1]*m[2][2];
  result -= m[0][0]*m[1][2]*m[2][1];
  result -= m[0][1]*m[1][0]*m[2][2];
  result += m[0][1]*m[1][2]*m[2][0];
  result += m[0][2]*m[1][0]*m[2][1];
  result -= m[0][2]*m[1][1]*m[2][0];
  // there's probably a more elegant way to do that
  return result;
}

void SymmetryMatrices::showMatrix(const Matrix &m) const {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      cout << setw(8) << m[i][j];
    }
    cout << endl;
  }
  return;
}

Matrix SymmetryMatrices::invert(const Matrix &m) const {
  double det = determinant(m);
  if (det == 0) {
    cerr << "ERROR: Symmetry matrix has no inverse!\n";
    showMatrix(m);
    exit(EXIT_FAILURE);
  }
  Matrix result;
  result.resize(3);
  for (int i = 0; i < 3; i++) {
    result[i].resize(3);
  }
  result[0][0] = m[1][1]*m[2][2] - m[1][2]*m[2][1];
  result[0][1] = m[0][2]*m[2][1] - m[0][1]*m[2][2];
  result[0][2] = m[0][1]*m[1][2] - m[0][2]*m[1][1];
  result[1][0] = m[1][2]*m[2][0] - m[1][0]*m[2][2];
  result[1][1] = m[0][0]*m[2][2] - m[0][2]*m[2][0];
  result[1][2] = m[0][2]*m[1][0] - m[0][0]*m[1][2];
  result[2][0] = m[1][0]*m[2][1] - m[1][1]*m[2][0];
  result[2][1] = m[0][1]*m[2][0] - m[0][0]*m[2][1];
  result[2][2] = m[0][0]*m[1][1] - m[0][1]*m[1][0];
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      result[i][j] /= det;
    }
  }
  return result;
}

// TODO: It may make sense to put together a lookup table of the 
// similarity-transformed matrices, products, etc.  I'm not sure how
// much of a speedup this might provide. 

Matrix SymmetryMatrices::globalRotationTransform(const Matrix &m) const {
  Matrix result = mult(globalRotation,mult(m,invert(globalRotation)));
  return result;
}

const Matrix SymmetryMatrices::getMatrix(int n) const {
  return globalRotationTransform(matrices->at(n));
}

const Matrix SymmetryMatrices::getInverse(int n) const {
  return globalRotationTransform(inverses->at(n));
}

const Matrix SymmetryMatrices::getProduct(int n, int m) const {
  return globalRotationTransform(products->at(n)[m]);
}

const Vec3 SymmetryMatrices::getOffsetPair(int n, int m) const {
  // TODO: My global rotation correction might do horrible things 
  // to symmetry with offset vectors; this definitely should be
  // checked.
  return offsetPairs->at(n)[m];
}

size_t SymmetryMatrices::size() const {
  return matrices->size();
}


