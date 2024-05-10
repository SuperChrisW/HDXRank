#include "AmberTrajectoryInput.h"
#include <cstdlib>
#include <sstream>
#include <iomanip>
using namespace std;

AmberTrajectoryInput::AmberTrajectoryInput() :
  frame(0),
  verbose(true)
{
}

AmberTrajectoryInput::AmberTrajectoryInput(
    const std::string &filename,
    const AmberPrmtop &prmtop ) :
    verbose( true ) {
  open( filename, prmtop );
}
AmberTrajectoryInput::AmberTrajectoryInput(
    const std::string &filename,
    int natoms_,
    bool ifbox_ ) :
    verbose( true ) {
  open( filename, natoms_, ifbox_ );
}

AmberTrajectoryInput::~AmberTrajectoryInput()
{
}

void AmberTrajectoryInput::rewind() {
  trajfile.seekg( 0, ios::beg );
  if ( !getline( trajfile, titlestring ) ) {
    cout << "Error: could not read trajectory title line" << endl;
  }
  frame = 0;
}

void AmberTrajectoryInput::open( const std::string &filename, const AmberPrmtop &prmtop ) {
  nAtoms = prmtop.natom;
  ifbox = prmtop.ifbox;
  open( filename );
}

void AmberTrajectoryInput::open( const std::string &filename, int natoms_, bool ifbox_ ) {
  nAtoms = natoms_;
  ifbox = ifbox_;
  open ( filename );
}

void AmberTrajectoryInput::open( const std::string &filename ) {
  trajfile.open( filename.c_str(), ios::in );
  if ( !trajfile ) {
    cout << "Error: could not open " << filename << endl;
    exit(0);
  }
  if ( !getline( trajfile, titlestring ) ) {
    cout << "Error: could not read trajectory title line from " << filename << endl;
    exit(0);
  }
  frame = 0;
}

void AmberTrajectoryInput::close() {
  trajfile.close();
  trajfile.clear();
}

bool AmberTrajectoryInput::readNextFrame( vector<Vec3> &coords ) {
  bool good;
  readNextFrame( coords, good );
  return good;
}

void AmberTrajectoryInput::readNextFrame( vector<Vec3> &coords, bool &good ) {
  double box0, box1, box2; //these are throw-away variables
  readNextFrame( coords, box0, box1, box2, good );

  if ( verbose ) {
    if ( good ) {
      if ( (frame%50)==0 ) {
        cout << "\n" << setw(6) << right << frame << flush;
      }
      cout << "." << flush;
    }
    else {
      cout << "\n" << frame << " frames successfully read from trajectory." << endl;
    }
  }
  frame++;
}

void AmberTrajectoryInput::readNextFrame( vector<Vec3> &coords, double &box0, double &box1, double &box2, bool &good ) {

  box0 = 0;
  box1 = 0;
  box2 = 0;
  string currentline;

  coords.resize( nAtoms );
  int atom = 0;
  int cartesianIndex = 0;
  int dataCount = 0;
  const int dataLimit = ifbox ? (nAtoms + 1)*3 : nAtoms*3;
  double tempdouble;
  while ( dataCount < dataLimit && trajfile >> tempdouble ) {

    //store the field
    if ( atom < nAtoms ) {
      switch ( cartesianIndex ) {
      case 0:
        coords[atom].x = tempdouble;
        break;
      case 1:
        coords[atom].y = tempdouble;
        break;
      case 2:
        coords[atom].z = tempdouble;
        break;
      }
    }
    else if ( ifbox ) {
      switch ( cartesianIndex ) {
      case 0:
        box0 = tempdouble;
        break;
      case 1:
        box1 = tempdouble;
        break;
      case 2:
        box2 = tempdouble;
        break;
      }
    }

    //increment counters
    dataCount++;
    cartesianIndex++;
    if ( cartesianIndex == 3 ) {
      cartesianIndex = 0;
      atom++;
    }
  }

  good = ( dataCount == dataLimit );
}
