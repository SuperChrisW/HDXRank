// Author name: Craig Jolley
// Created:     19 Feb 2008

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include "EZDMap.h"


///////////////////////////////////////////////////////////////////////////////
// Description:
//   Constructor for EZDMap class.  Loads in the information given in the file 
//   header (unit cell parameters, origin, extent, number of grid points, and 
//   scaling factor), loads the ED information from the map, calculates the 
//   unit cell vectors and the real-space offset from the origin.
// Parameters:
//   char *fileName -- name of the EZD-format file containing the ED map
///////////////////////////////////////////////////////////////////////////////
EZDMap::EZDMap(std::string fileName) {
  using namespace std;
  char err1[36] = "Error in file header: problem with ";
  char err2[22] = " line.\n";
  ifstream inFile;
  inFile.open(fileName.c_str());
  // First, make sure the file is kosher
  if (!inFile.is_open()) {
    cerr << "Failed to open file " << fileName << " for reading.\n";
    exit(EXIT_FAILURE);
  }
  string test;
  inFile >> test;  // first line should be "EZD_MAP"
  if (test!="EZD_MAP") {
    cerr << "Error in file header: not an EZD map!\n";
    inFile.close();
    exit(EXIT_FAILURE);
  }
  // Extract information from the file header
  cout << "Loading header information from " << fileName << " ...\n";
  inFile.get();                // throw away \n at end of first line
  inFile.getline(comment,80);  // second line is a comment
  inFile >> test;              // third line should begin with "CELL"
  if (test!="CELL") {
    cerr << err1 << "CELL" << err2;
    inFile.close();
    exit(EXIT_FAILURE);
  }
  inFile >> a;                // get unit cell parameters
  inFile >> b;
  inFile >> c;
  inFile >> alpha;
  inFile >> beta;
  inFile >> gamma;
  inFile >> test;              // fourth line should begin with "ORIGIN"
  if (test != "ORIGIN")	{
    cerr << err1 << "ORIGIN" << err2;
    inFile.close();
    exit(EXIT_FAILURE);
  }
  inFile >> origin.i;         // get origin
  inFile >> origin.j;
  inFile >> origin.k;
  inFile >> test;             // fifth line should begin with "EXTENT"
  if (test != "EXTENT") {
    cerr << err1 << "EXTENT" << err2;
    inFile.close();
    exit(EXIT_FAILURE);
  }
  inFile >> extent.i;        // get grid extent
  inFile >> extent.j;
  inFile >> extent.k;
  inFile >> test;             // sixth line should begin with "GRID"
  if (test != "GRID")	{
    cerr << err1 << "GRID" << err2;
    inFile.close();
    exit(EXIT_FAILURE);
  }
  inFile >> num.i;        // get number of grid points in each direction
  inFile >> num.j;
  inFile >> num.k;
  inFile >> test;             // seventh line should begin with "SCALE"
  if (test != "SCALE") {
    cerr << err1 << "SCALE" << err2;
    inFile.close();
    exit(EXIT_FAILURE);
  }
  inFile >> scale;      // get scaling factor
  // Load in grid data from file
  cout << "Loading grid data from " << fileName << ":   0%";
  dataSize = extent.i*extent.j*extent.k;  // # of points on grid
  gridData = new double[dataSize];
  inFile >> test;
  if (test!="MAP") {  // eighth line should begin with "MAP"
    cerr << err1 << "MAP" << err2;
    inFile.close();
    exit(EXIT_FAILURE);
  }
  int lastPC = 0; // previous percentage complete; see below
  for (long n = 0; n < dataSize; n++) {
    double inputNumber;
    inFile >> inputNumber;
    gridData[n] = inputNumber / scale;
    // Progress report for loading big files
    if (n % 1000 == 0) {
      int percentComplete = int (100.0 * n / dataSize);
      if ((percentComplete % 5 == 0)&&(percentComplete > lastPC)) {
        cout << "\b\b\b\b" << setw(3) << percentComplete << '%' << flush;
	lastPC = percentComplete;
      }
    }
  }
  cout << endl;
  // Calculate unit-cell vectors i, j, and k
  double iScale = a * extent.i / num.i;
  double jScale = b * extent.j / num.j;
  double kScale = c * extent.k / num.k;
  // for orthorhombic special case, avoid round-off errors
  if (alpha == 90.0 && beta == 90.0 && gamma == 90.0) {
    i = Vec3(iScale,0.0,0.0);
    j = Vec3(0.0,jScale,0.0);
    k = Vec3(0.0,0.0,kScale);
  } else {
    double twoPi360 = 0.0174532925199;
    double radAlpha, radBeta, radGamma;
    radAlpha = twoPi360 * alpha;
    radBeta = twoPi360 * beta;
    radGamma = twoPi360 * gamma;
    i = Vec3(iScale, 0.0, 0.0);
    j = Vec3(jScale*cos(radGamma), jScale*sin(radGamma), 0.0);
    k.x = kScale*cos(radBeta);
    k.y = kScale*(cos(radAlpha) - cos(radBeta)*cos(radGamma))/sin(radGamma);
    k.z = sqrt(pow(kScale,2) - pow(k.x,2) - pow(k.y,2));
  }
  // Now calculate the real-space offset
  offset.x = (i.x * (origin.i) + j.x * (origin.j) + k.x * origin.k)/extent.i;
  offset.y = (i.y * (origin.i) + j.y * (origin.j) + k.y * origin.k)/extent.j;
  offset.z = (i.z * (origin.i) + j.z * (origin.j) + k.z * origin.k)/extent.k;
  displayHeader();
  inFile.close();
  return;
}
///////////////////////////////////////////////////////////////////////////////
// Description:
//   Displays the information contained in the EZD file header
///////////////////////////////////////////////////////////////////////////////
void EZDMap::displayHeader() {
  using namespace std;
  cout << "  Information from EZD header file:" << endl << endl;
  cout << comment << endl;
  cout << "  a = " << a << ", b = " << b << ", c = " << c << endl;
  cout << "  alpha = " << alpha << ", beta = " << beta << ", gamma = " << gamma << endl;
  cout << "  Origin: " << origin.i << ", " << origin.j << ", " << origin.k << endl;
  cout << "  Grid extent is " << extent.i << " x " << extent.j << " x " << extent.k << endl;
  cout << "  Num of grid points: " << num.i << " x " << num.j << " x " << num.k << endl;
  cout << "  Scaling factor = " << scale << "\n\n";
  return;
}



