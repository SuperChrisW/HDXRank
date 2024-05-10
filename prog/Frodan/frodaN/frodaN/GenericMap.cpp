// Author name: Craig Jolley
// Created:     June 2006

// NOTE: This code was originally written (in 2006) to work with classic
//       FRODA and was in Froda_ED.cpp.  I've started (Feb 2008) to adapt 
//       it for use with newsim.  Not all of the functionality in the FRODA 
//       version will be included for now; I'm omitting the TheoMap and 
//       TrimMap classes.  Consequently, if you want to generate a 
//       theoretical cryo-EM map or trim an existing map, use 
//       classic FRODA for those steps and do the fitting in newsim.

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <vector>
#include <ctime>
#include "tricubic.h"
#include "GenericMap.h"

////////////////////////////////////////////////////////////////////////////////
// Description:  
//   Constructor for GenericMap class.  
////////////////////////////////////////////////////////////////////////////////
GenericMap::GenericMap(){
  diffMap = NULL;  // instantiate if needed
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Destructor for GenericMap class.
////////////////////////////////////////////////////////////////////////////////
GenericMap::~GenericMap(){
  if (diffMap != NULL) { delete [] diffMap; }
}


///////////////////////////////////////////////////////////////////////////////
// Description:
//   Given crystallographic grid indices i,j,k, the location of the 
//   corresponding point in the 1-D array GenericMap::gridData is given by:
//   n = i + extent.i*j + extent.i*extent.j*k
//   This function does the reverse transformation
// Parameters:
//   int n -- a position in the 1-D array GenericMap::gridData
// Return Value List:
//   Returns a gridPoint structure containing the i,j,k indices corresponding 
//   to the 1-D index n
///////////////////////////////////////////////////////////////////////////////
gridPoint GenericMap::gridIndex(int n) {
  gridPoint temp;
  temp.i = n % extent.i;
  temp.j = (n % (extent.i*extent.j) - temp.i)/extent.i;
  temp.k = (n - temp.j * extent.i - temp.i)/(extent.i*extent.j);
  return temp;
}

///////////////////////////////////////////////////////////////////////////////
// Description:
//   Given a grid point i,j,k, finds its location in Cartesian coordinates.
// Parameters:
//   gridPoint gp -- the grid point whose position is sought
// Return Value list:
//   Returns a Vector object containing the location of gp in Cartesian 
//   coordinates.
///////////////////////////////////////////////////////////////////////////////
Vec3 GenericMap::gridPointToVector(gridPoint gp) {
  Vec3 temp(0.0, 0.0, 0.0);
  temp.x += gp.i * i.x / extent.i;
  temp.x += gp.j * j.x / extent.j;
  temp.x += gp.k * k.x / extent.k;
  temp.y += gp.i * i.y / extent.i;
  temp.y += gp.j * j.y / extent.j;
  temp.y += gp.k * k.y / extent.k;
  temp.z += gp.i * i.z / extent.i;
  temp.z += gp.j * j.z / extent.j;
  temp.z += gp.k * k.z / extent.k;
  temp += offset;
  return temp;
}


///////////////////////////////////////////////////////////////////////////////
// Description:
//   After the map has been initialized, call this to output the map in EZD
//   format.
// Parameters:
//   string fileName -- the filename to output
///////////////////////////////////////////////////////////////////////////////
void GenericMap::writeEZD(std::string fileName) {
  using namespace std;
  ofstream outFile;
  outFile.open(fileName.c_str());
  if (!outFile.is_open()) {
    cerr << "Failed to output EZD map to " << fileName << endl;
    exit(EXIT_FAILURE);
  }
  cout << "Writing EZD map to " << fileName << endl;
  outFile << "EZD_MAP\n";
  outFile << "! Generated from theoretical electron density map.\n";
  outFile << "CELL " << a << ' ' << b << ' ' << c << ' ';
  outFile << alpha << ' ' << beta << ' ' << gamma << endl;
  outFile << "ORIGIN " << origin.i << ' ' << origin.j << ' ' << origin.k << endl;
  outFile << "EXTENT " << extent.i << ' ' << extent.j << ' ' << extent.k << endl;
  outFile << "GRID " << extent.i << ' ' << extent.j << ' ' << extent.k << endl;
  // find maximum value of electron density
  int gridSize = extent.i*extent.j*extent.k;
  double maxDensity = 0;
  for (long n = 0; n < gridSize; n++) {
    if (gridData[n] > maxDensity) {
      maxDensity = gridData[n];
    }
  }
  double scale = 9999.0 / maxDensity;
  outFile << "SCALE ";
  outFile << scientific << setprecision(8) << uppercase << scale << endl;
  outFile << "MAP\n";
  for (long n = 0; n < gridSize; n++) {
    outFile << fixed << setprecision(1) << gridData[n] * scale;
    if (((n+1) % 7 == 0)||n==gridSize-1) {  // seven numbers per line
      outFile << endl;
    }
    else {
      outFile << ' ';
    }
  }
  outFile << "END\n";
  outFile.close();
  return;
}

///////////////////////////////////////////////////////////////////////////////
// Description:
//   Returns the determinant of the 3x3 matrix whose rows are formed by v1,
//   v2, and v3.  
///////////////////////////////////////////////////////////////////////////////
double GenericMap::det3(Vec3 v1, Vec3 v2, Vec3 v3) {
  double result = 0.0;
  result += v1.x * (v2.y * v3.z - v2.z * v3.y);
  result -= v1.y * (v2.x * v3.z - v2.z * v3.x);
  result += v1.z * (v2.x * v3.y - v2.y * v3.x);
  return result;
}

///////////////////////////////////////////////////////////////////////////////
// Description:
//   Returns the "grid point" corresponding to a real space vector.  The real
//   grid points have integer indices; in general this will correspond to a 
//   point that is intermediate between eight defined grid points.  The x,y,z
//   components of the return value correspond to the i,j,k indices of the 
//   grid point.  
//   NOTE: It might seem wasteful to make these calls to det3() every time, 
//   since they're doing a lot of the same operations over and over again.
//   I tried calculating these at the beginning of a run and storing the 
//   results for later use; there was no noticeable speedup.
//////////////////////////////////////////////////////////////////////////////
Vec3 GenericMap::vectorToGridPoint(Vec3 v) {
  using namespace std;
  Vec3 result;
  v -= offset;
  double det_ijk = det3(i,j,k);
  result.x = extent.i * det3(v,j,k) / det_ijk;
  result.y = extent.j * det3(v,k,i) / det_ijk;
  result.z = extent.k * det3(v,i,j)  / det_ijk; 
  return result;
}

///////////////////////////////////////////////////////////////////////////////
// Description:
//   Basis polynomials for cubic interpolation and their derivatives.  
//   poly1(0,q) = 1 if q = -1 and 0 at q =  0, 1, 2
//   poly2(1,q) = 1 if q = 0  and 0 at q = -1, 1, 2
//   poly3(2,q) = 1 if q = 1  and 0 at q = -1, 0, 2
//   poly4(3,q) = 1 if q = 2  and 0 at q = -1, 0, 1
//   This means that a linear combination of these four polynomials will be 
//   the unique cubic polynomial determined by its values at {-1,0,1,2}.
//   Not used in the new cubic interpolation scheme.
//////////////////////////////////////////////////////////////////////////////
/*double GenericMap::poly(int n, double q) { 
  switch(n) {
    case 0 : return -q*(q-1)*(q-2)/6.0; 
    case 1 : return 3*(q+1)*(q-1)*(q-2)/6.0;
    case 2 : return -3*(q+1)*q*(q-2)/6.0;
    case 3 : return (q+1)*q*(q-1)/6.0;
  }
  return 0;
}
double GenericMap::dpoly(int n, double q) { 
  switch(n) {
    case 0 : return (-3*q*q + 6*q - 2)/6.0;
    case 1 : return (9*q*q - 12*q - 3)/6.0;
    case 2 : return (-9*q*q + 6*q + 6)/6.0; 
    case 3 : return (3*q*q-1)/6.0;
  }
  return 0;
  }*/

///////////////////////////////////////////////////////////////////////////////
// Description:
//   Given a real-space vector, returns the gradient of the map at that point.
///////////////////////////////////////////////////////////////////////////////
Vec3 GenericMap::mapGradient(const Vec3 &pos) {
  using namespace std;
  if (alpha != 90.0 || beta != 90.0 || gamma != 90.0) {
    std::cerr << "ERROR: Gradient of map can (currently) only be calculated for " 
	 << "orthorhombic unit cells (alpha = beta = gamma = 90 degrees).  "
	 << "You'll need to revise GenericMap::gradient() to account for a non-"
	 << "orthorhombic unit cell.  Lucky you.\n";
    exit(EXIT_FAILURE);
  }
  Vec3 gridLocation = vectorToGridPoint(pos);
  // the interpolation is over the region bounded by the points 
  // 000,100,010,001,011,101,110,111.  The variables defied below are the
  // coordinates of the 000 point
  gridPoint o;
  o.i = (int) gridLocation.x;
  o.j = (int) gridLocation.y;
  o.k = (int) gridLocation.z;
  //cout << "gridLocation: " << gridLocation << endl;
  //cout << "origin: " << origin.i << ',' << origin.j << ',' << origin.k << endl;
  // In addition to the eight grid points mentioned above, the shell formed by
  // the nearest neighbors of these grid points (with indices ranging from -1 
  // to 2) is used to provide 64 constraints for the tricubic interpolation.  
  // C++ likes to have arrays start with 0, so the indices I actually use will 
  // run from 0 to 3.
  double local[4][4][4]; 
  for (int l = 0; l < 4; l++) {
    for (int m = 0; m < 4; m++) {
      for (int n = 0; n < 4; n++) {
        if (o.i+l-1 < 0 || o.i+l-1 > extent.i ||
            o.j+m-1 < 0 || o.j+m-1 > extent.j ||
            o.k+n-1 < 0 || o.k+n-1 > extent.k) {
          local[l][m][n] = 0;
	} else {
          local[l][m][n] = gridData[(o.i+l-1) + extent.i*(o.j+m-1) + 
            extent.i*extent.j*(o.k+n-1)];
	}
      }
    }
  }
  Vec3 originPos = gridPointToVector(o);
  Vec3 scaledPos;
  // The following lines assume that the unit cell is orthorhombic -- this
  // means that i points solely along the x direction, j along the y, and 
  // k along the z direction.
  scaledPos.x = (pos.x - originPos.x) / (i.x / extent.i);
  scaledPos.y = (pos.y - originPos.y) / (j.y / extent.j);
  scaledPos.z = (pos.z - originPos.z) / (k.z / extent.k);
  // scaledPos is now in terms of a set of scaled coordinates -- its components take values on the
  // interval [0,1)
  Vec3 result(0,0,0);
  /*
  // This is my ghetto cubic interpolation; the C1 interpolation may work better
  for (int l = 0; l < 4; l++) {
    for (int m = 0; m < 4; m++) {
      for (int n = 0; n < 4; n++) {
        // It looks like calculating these three polynomials once instead of 
        // twice gets a small speedup (about 8% in my test case)
        double poly_lx = poly(l,scaledPos.x);
        double poly_my = poly(m,scaledPos.y);
        double poly_nz = poly(n,scaledPos.z);
        result.x += dpoly(l,scaledPos.x) * poly_my * poly_nz * local[l][m][n];
        result.y += poly_lx * dpoly(m,scaledPos.y) * poly_nz * local[l][m][n];
        result.z += poly_lx * poly_my * dpoly(n,scaledPos.z) * local[l][m][n];
      }
    }
  }
  // now re-scale the gradient the same way coordinates were re-scaled
  result.x /= (i.x / extent.i);
  result.y /= (j.y / extent.j);  
  result.z /= (k.z / extent.k);
  */

  // use the C1 interpolation scheme of Lekien et al
  // set up arrays
  double f[8];
  double dfdx[8];     double dfdy[8];     double dfdz[8];
  double d2fdxdy[8];  double d2fdxdz[8];  double d2fdydz[8];
  double d3fdxdydz[8];
  int index = 0;
  // this calculation will be the same for each iteration; it may be worthwhile to 
  // calculate and store the results
  for (int z = 1; z <= 2; z++) {
    for (int y = 1; y <= 2; y++) {
      for (int x = 1; x <= 2; x++) {
        f[index] = local[x][y][z];
        dfdx[index] = 0.5*(local[x+1][y][z] - local[x-1][y][z])*(i.x / extent.i);
	dfdy[index] = 0.5*(local[x][y+1][z] - local[x][y-1][z])*(j.y / extent.j);
        dfdz[index] = 0.5*(local[x][y][z+1] - local[x][y][z-1])*(k.z / extent.k);
        d2fdxdy[index] = 0.25*(local[x-1][y-1][z] - local[x+1][y-1][z] - local[x-1][y+1][z] 
			       + local[x+1][y+1][z])*(i.x / extent.i)*(j.y / extent.j);
        d2fdxdz[index] = 0.25*(local[x-1][y][z-1] - local[x+1][y][z-1] - local[x-1][y][z+1] 
                               + local[x+1][y][z+1])*(i.x / extent.i)*(k.z / extent.k);
        d2fdydz[index] = 0.25*(local[x][y-1][z-1] - local[x][y-1][z+1] - local[x][y+1][z-1] 
                               + local[x][y+1][z+1])*(j.y / extent.j)*(k.z / extent.k);
        d3fdxdydz[index] = 0.125*(local[x+1][y+1][z+1] + local[x+1][y-1][z-1] + local[x-1][y+1][z-1] 
                                + local[x-1][y-1][z+1] - local[x-1][y+1][z+1] - local[x+1][y-1][z+1]
				- local[x+1][y+1][z-1] - local[x-1][y-1][z-1])
	                        *(i.x / extent.i)*(j.y / extent.j)*(k.z / extent.k);
                                    
        index++;
        // Multiplying by the grid spacings here might seem a little unusual.  This is needed because
        // just taking sums and differences of the map values would give you approximations to  
        // derivatives in the *unscaled* coordinate system, while the tricubic interpolator needs
        // numbers in the *scaled* coordinate system.  If the unscaled spatial variable is x and
        // the scaled one is x' = x/s, where s is the grid spacing, then the scaled derivative
        // will be df/dx' = df/dx * dx/dx' = df/dx * s
      }
    }
  }
  // call tricubic library functions
  double a[64];
  tricubic_get_coeff(a,f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz);
  // this step will differ for each iteration
  result.x = tricubic_eval(a,scaledPos.x,scaledPos.y,scaledPos.z,1,0,0) / (i.x / extent.i);
  result.y = tricubic_eval(a,scaledPos.x,scaledPos.y,scaledPos.z,0,1,0) / (j.y / extent.j);
  result.z = tricubic_eval(a,scaledPos.x,scaledPos.y,scaledPos.z,0,0,1) / (k.z / extent.k);
  // This is exactly the reverse of the situation we had above; we want to return the gradient in
  // the *unscaled* coordinate system, so we need df/dx = df/dx' * dx'/dx = df/dx' * 1/s
  return result;
}
  
///////////////////////////////////////////////////////////////////////////////
// Description:
//   Calculates the correlation function C between the an electron density map
//   and an EMStructure object.  
// Parameters:
//   EMStructure &pdb -- the EMStructure object whose correlation with the EM
//                       map is to be computed
// Return Value List:
//   The correlation is a real number between 0.0 and 1.0, with a score of 1.0 
//   indicating perfect correlation and 0.0 indicating no overlap at all 
//   between the two map densities. 
///////////////////////////////////////////////////////////////////////////////
double GenericMap::correlate(EMStructure &em, double &denom1, double &denom2) {
  using namespace std;
  double cNumerator = 0.0;     // numerator of C
  double d1 = 0.0;       // sums in denominator of C
  double d2 = 0.0;
  em.updateGrid();
  for (long n = 0; n < dataSize; n++) {
    double localDensity = em.densityAt(gridPointToVector(gridIndex(n)));
    cNumerator += gridData[n]*localDensity;
    d1 += gridData[n]*gridData[n];
    // It might seem like you could save time by saving the value of
    // cDenominator1 because it's conformation-invariant, but it actually
    // doesn't make a difference in speed.  You're welcome.
    d2 += localDensity*localDensity;
  }  
  double result = cNumerator/sqrt(d1*d2);
  denom1 = d1;
  denom2 = d2;
  return result;
}

///////////////////////////////////////////////////////////////////////////////
// Description:
//    Calculates a scaled difference map needed by cGradient().
// Parameters:
//    EMStructure &em -- EMStructure to be compared to map
///////////////////////////////////////////////////////////////////////////////
void GenericMap::calcDiffMap(EMStructure &em) {
  if (diffMap == NULL) { diffMap = new double[dataSize]; }
  double denom1, denom2;
  double corr = correlate(em,denom1,denom2);
  double sd1 = sqrt(denom1);
  double sd2 = sqrt(denom2);
  for (long n = 0; n < dataSize; n++) {
    double localDensity = em.densityAt(gridPointToVector(gridIndex(n)));
    diffMap[n] = gridData[n]/sd1 - corr*localDensity/sd2;
  } 
  std::cout << "  Correlation: " << corr << std::endl;
  return;
}

///////////////////////////////////////////////////////////////////////////////
// Description:
//   Returns the indices of points in the grid near a specified point in 
//   space.  Needed by cGradient().
// Parameters:
//   long n -- index (in gridData[]) of the point in question
//   double cutoff -- distance within which all points should be included
// Return value:
//   STL vector containing a list of all nearby points, as gridPoint objects
///////////////////////////////////////////////////////////////////////////////
std::vector <gridPoint> GenericMap::nearbyPoints(Vec3 v, double cutoff) {
  using namespace std;
  // the following assumes a cubic grid
  Vec3 spacing(i.x/extent.i, j.y/extent.j, k.z/extent.k);
  gridPoint max;
  max.i = (int) ceil(cutoff/spacing.x)+1;
  max.j = (int) ceil(cutoff/spacing.y)+1;
  max.k = (int) ceil(cutoff/spacing.z)+1;
  Vec3 nearGrid = vectorToGridPoint(v);
  gridPoint gp;
  gp.i = (int) (nearGrid.x + 0.5);
  gp.j = (int) (nearGrid.y + 0.5);
  gp.k = (int) (nearGrid.z + 0.5);
  gridPoint start, finish;
  if (gp.i <= max.i) { start.i = 0; } else { start.i = gp.i - max.i; }
  if (gp.j <= max.j) { start.j = 0; } else { start.j = gp.j - max.j; }
  if (gp.k <= max.k) { start.k = 0; } else { start.k = gp.k - max.k; }
  if (gp.i >= extent.i - max.i) { finish.i = extent.i; } else { finish.i = gp.i + max.i; }
  if (gp.j >= extent.j - max.j) { finish.j = extent.j; } else { finish.j = gp.j + max.j; }
  if (gp.k >= extent.k - max.k) { finish.k = extent.k; } else { finish.k = gp.k + max.k; }
  gridPoint iter;
  std::vector <gridPoint> result;
  for (iter.i = start.i; iter.i < finish.i; iter.i++) {
    for (iter.j = start.j; iter.j < finish.j; iter.j++) {
      for (iter.k = start.k; iter.k < finish.k; iter.k++) {
        if (v.dist2(gridPointToVector(iter)) < cutoff*cutoff) {
          result.push_back(iter);
        }
      }
    }
  }
  return result;
}

///////////////////////////////////////////////////////////////////////////////
// Description:
//   Calculates the gradient of the correlation function with respect to the 
//   position of a particular atom -- this will provide the force applied
//   during correlation-based perturbation
// Parameters:
//   EMStructure &em -- the EMStructure object whose correlation with the EM
//                       map is to be computed
//   int atomID -- the ID (in rigidUnitSystem) of the atom in question
//   double corr -- the correlation needs to be calculated first and is rather
//                  computationally expensive; it's best to find it once and
//                  then use it repeatedly for subsequent calls of cGradient
// Return Value List:
//   Returns a real-space vector with the gradient.
///////////////////////////////////////////////////////////////////////////////
Vec3 GenericMap::cGradient(EMStructure &em, int atomID) {
  Vec3 result(0,0,0);
  Vec3 atomVec = em.atomPos(atomID);
  std::vector <gridPoint> np = nearbyPoints(atomVec,em.cutoff);
  // I also tried a version of nearbyPoint() that passes np by reference; it 
  // didn't improve performance.
  for (size_t i = 0; i < np.size(); i++) {
    Vec3 addMe = gridPointToVector(np[i]) - atomVec;
    addMe *= em.scaledDensity(addMe.norm2(),atomID) * diffMap[gridPointToArray(np[i])];
    result += addMe;
  } 
  return result;
}
  
  
