// Author name:  Craig Jolley
// Created:      19 Jun 2006

// NOTE: This code was originally written (in 2006) to work with classic
//       FRODA; I've started (Feb 2008) to adapt it for use with newsim.
//       Not all of the functionality in the FRODA version will be included
//       for now; I'm omitting the TheoMap and TrimMap classes.
//       Consequently, if you want to generate a theoretical cryo-EM map or 
//       trim down an existing map, use classic FRODA for those
//       steps and do the fitting in newsim.

#ifndef GENERIC_MAP_H_
#define GENERIC_MAP_H_

#include <string>
#include <vector>
#include "Vec3.h"
#include "EMStructure.h"

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Indices of a single point in the electron density grid.  Used internally
//   by the GenericMap class.
////////////////////////////////////////////////////////////////////////////////
struct gridPoint {  // single point on grid
  int i;
  int j;
  int k;
};

////////////////////////////////////////////////////////////////////////////////
// Description:
//   An abstract electron density map class.  This contains the features that
//   should be common to all types of electron density maps, e.g. storing the
//   electron density data, converting between its crystallographic grid 
//   coordinates and 3D Cartesian coordinates, and calculating its correlation
//   with an EMStructure object.
////////////////////////////////////////////////////////////////////////////////
class GenericMap {
protected:
  double a; // unit cell parameters
  double b;
  double c;
  double alpha;
  double beta;
  double gamma;
  int dataSize;      // equal to extent.i*extent.j*extent.k
  Vec3 i, j, k;    // unit cell vectors
  Vec3 offset;     // real-space offset from origin
  gridPoint origin;  // grid origin
  double *gridData;  // array containing map data
  double *diffMap; // difference map for cc gradient perturbation
  gridPoint extent;   // grid extent
private:
  Vec3 vectorToGridPoint(Vec3 v);
  long gridPointToArray(gridPoint gp) {
    return gp.i + extent.i*gp.j + extent.i*extent.j*gp.k; }
  Vec3 gridPointToVector(gridPoint gp);  // convert i,j,k to x,y,z
  //double poly(int n,double q);    // basis polynomials for cubic interpolation
  //double dpoly(int n, double q);  // and their derivatives
  double correlate(EMStructure &em, double &denom1, double &denom2);
  // calculate real-space correlation; returns the two denominator terms,
  // since they're also needed by cGradient
  std::vector <gridPoint> nearbyPoints(Vec3 v, double cutoff);
  //returns points within cutoff of point in space
public:
  GenericMap();
  ~GenericMap();
  gridPoint getExtent() {return extent;}
  void writeEZD(std::string fileName); // output EZD map
  Vec3 mapGradient(const Vec3 &pos);
  void calcDiffMap(EMStructure &em); 
  // calculates the difference map needed by cGradient()
  Vec3 cGradient(EMStructure &em, int atomID); 
  // gradient of correlation with respect to an atom position
protected:
  double det3(Vec3 v1, Vec3 v2, Vec3 v3); // determinant of 3x3 matrix
  gridPoint gridIndex(int n);
   // converts a position in the 1D data array to a position on the 3D grid
};

#endif
