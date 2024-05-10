#include "AmberTrajectory.h"
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <iomanip>

using namespace std;

AmberTrajectory::AmberTrajectory() :
  readyToAppend(false)
{
}

AmberTrajectory::~AmberTrajectory()
{
}

void AmberTrajectory::initializeOutputFile( const string &filename_, bool box ) {
  filename = filename_;
  ofstream outfile( filename.c_str(), ios::out );
  outfile << '\n';
  outfile.flush();
  outfile.close();
  ifbox = box;
  width = 8;
  precision = 3;
  nCols = 10;
  readyToAppend = true;
}

void AmberTrajectory::writeRestartFile( const string &filename_, bool box, const vector<Vec3> &coords ) {
  filename = filename_;
  ofstream outfile( filename.c_str(), ios::out );
  outfile << '\n';
  outfile << setw(6) << coords.size() << '\n';
  outfile.flush();
  outfile.close();
  ifbox = box;
  width = 12;
  precision = 7;
  nCols = 6;
  readyToAppend = true;
  append( coords );
  readyToAppend = false;
}

void AmberTrajectory::append( const vector<Vec3> &coords ) {
  if ( !readyToAppend ) {
    cout << "AmberTrajectory: Cannot append data.  File is not yet initialized." << endl;
    exit(0);
  }

  ofstream outfile( filename.c_str(), ios::app );
  outfile << showpoint << setiosflags(ios::fixed) << setprecision( precision );
  int j = 0;
  for ( size_t i = 0; i < coords.size(); i++ ) {
    outfile << setw( width ) << coords[i].x;
    if ( ++j == nCols ) { j = 0; outfile << '\n'; }
    outfile << setw( width ) << coords[i].y;
    if ( ++j == nCols ) { j = 0; outfile << '\n'; }
    outfile << setw( width ) << coords[i].z;
    if ( ++j == nCols ) { j = 0; outfile << '\n'; }
  }

  //If AMBER prmtop requires a box, write out dummy box
  //information (zeros)
  if ( ifbox ) {
    outfile << setw( width ) << 0.0;
    if ( ++j == nCols ) { j = 0; outfile << '\n'; }
    outfile << setw( width ) << 0.0;
    if ( ++j == nCols ) { j = 0; outfile << '\n'; }
    outfile << setw( width ) << 0.0;
    if ( ++j == nCols ) { j = 0; outfile << '\n'; }
  }

  // To properly separate the frames in the AMBER
  // trajectory file,
  // after writing all the coordinates, we must make sure
  // that the frame is ended with a newline.  If j==0,
  // then the required newline character has already
  // just been written.  If not, it will be added here.
  if ( j ) { j = 0; outfile << '\n'; }

  outfile.flush();
  outfile.close();
}
