#include "AmberRestartFileInput.h"
#include <sstream>
#include <cstdlib>
using namespace std;

AmberRestartFileInput::AmberRestartFileInput()
{
}

AmberRestartFileInput::~AmberRestartFileInput()
{
}

void AmberRestartFileInput::read( const std::string &filename, const AmberPrmtop &prmtop, std::vector<Vec3> &coords, bool &good ) {
  std::ifstream infile;
  infile.open( filename.c_str(), ios::in );
  if ( infile.bad() ) {
    cout << "Error in AmberRestartFileInput: could not open " << filename << endl;
    exit(0);
  }
  
  string currentline;

  //read first line, which is the title line
  getline( infile, currentline );

  good = infile.good();
  if ( !good ) {
    cout << "Error in AmberRestartFileInput: could not read title line " << filename << endl;
    exit(0);
  }

  double box0;
  double box1;
  double box2;
  
  //read second line, which contains the number of atoms
  getline( infile, currentline );
  int nAtoms; 
  std::stringstream ss;
  ss << currentline;
  ss >> nAtoms;
  good = nAtoms == prmtop.natom;
  if ( !good ) {
    cout << "Error in AmberRestartFileInput: nAtoms does not match prmtop nAtoms "  << endl;
    exit(0);
  }

  
  int ifbox = prmtop.ifbox;

  coords.resize( nAtoms );
  int atom = 0;
  int cartesianIndex = 0;
  int dataCount = 0;
  const int dataLimit = ifbox ? (nAtoms + 1)*3 : nAtoms*3;
  const int fieldWidth = 12;
  double tempdouble;
  while ( infile.good() && !infile.eof() && dataCount < dataLimit ) { //loop over lines
    getline( infile, currentline );
    int offset=0;
    while ( true ) { // loop over fields in the line
      //read the next field
      if ( offset + fieldWidth > static_cast<int> (currentline.size()) ) break;
      std::string fieldString = currentline.substr( offset, fieldWidth );
      if ( fieldString.find_first_not_of(" \n\t") == std::string::npos ) break;
      std::stringstream ss;
      ss << fieldString;
      ss >> tempdouble;
      
      //store the field
      if ( atom < nAtoms ) {
        if ( cartesianIndex == 0 ) coords[atom].x = tempdouble;
        if ( cartesianIndex == 1 ) coords[atom].y = tempdouble;
        if ( cartesianIndex == 2 ) coords[atom].z = tempdouble;
      }
      else if ( atom == nAtoms ) {
        if ( cartesianIndex == 0 ) box0 = tempdouble;
        if ( cartesianIndex == 1 ) box1 = tempdouble;
        if ( cartesianIndex == 2 ) box2 = tempdouble;        
      }
      else {
        good = false;
        cout << "Error: This should never happen" << endl;
        cout << "atom " << atom << " nAtoms " << nAtoms << endl;
        exit(0);
      }
      
      //increment counters
      dataCount++;
      cartesianIndex++;
      if ( cartesianIndex == 3 ) {
        cartesianIndex = 0;
        atom++;
      }
      offset+=fieldWidth;
    } 
  }
  
  if ( dataCount != dataLimit ) {
    good = false;
  }
}
