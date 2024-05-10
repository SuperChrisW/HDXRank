#include "TextFileInput.h"
#include "RigidUnitSystem.h"
//#include "Ropes.h"
#include "NeighborTable.h"
#include "PairTypeCutoffs.h"
#include "AmberPrmtop.h"
#include "TargetEnergy.h"
#include "TargetMap.h"
#include "HBManager.h"
#include "PHManager.h"
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
using namespace std;

TextFileInput::TextFileInput()
{
}

TextFileInput::~TextFileInput()
{
}

/*
Ropes *TextFileInput::buildRopes(
  string filename,
  const RigidUnitSystem *rigidUnitSystem )
{
  Ropes *ropes = new Ropes( rigidUnitSystem );

  ifstream pairfile( filename.c_str(), ios::in );
  if (!pairfile) {
    cout << "Could not open file: " << filename << endl;
    exit(1);
  }
  string currentline;
  int p1;
  int p2;
  double ropeLength;
  while (!pairfile.eof()) {
    getline( pairfile, currentline );
    if ( currentline.size() < 3 ) continue;
    stringstream ss;
    ss.str( currentline );
    ss >> p1;
    ss >> p2;
    ss >> ropeLength;
    ss.str("");
    ss.clear();
    if ( p2 < p1 ) swap( p1, p2 );
    if ( p1 == p2 ) continue;

    if ( rigidUnitSystem->doPointsBelongToSameRigidUnit( p1, p2 ) ) continue;
    //1-4 interactions are not explicitly excluded here, but it is assumed
    //that the list does not contain any
    ropes->insert(p1, p2, ropeLength);
  }
  return ropes;
}
*/

void TextFileInput::readConstraintsFile(
  std::string filename,
  vector<int>& firstatom,
  vector<int>& secondatom,
  vector<double>& dist ) {

  ifstream pairfile( filename.c_str(), ios::in );
  if (!pairfile) {
    cout << "Could not open file: " << filename << endl;
    exit(1);
  }
  string currentline;
  int p1;
  int p2;
  double cutoff;
  while (!pairfile.eof()) {
    getline( pairfile, currentline );
    if ( currentline.size() < 3 ) continue;
    stringstream ss;
    ss.str( currentline );
    ss >> p1;
    ss >> p2;
    ss >> cutoff;
    ss.str("");
    ss.clear();
    if ( p2 < p1 ) swap( p1, p2 );
    if ( p1 == p2 ) continue;

    firstatom.push_back( p1 );
    secondatom.push_back( p2 );
    dist.push_back( cutoff );
  }

}

void TextFileInput::readHBondFile_index0( string filename, HBManager *hbManager ) {
      cout << "Reading h-bonds from file " << filename << endl;
      cout << "  Interpreting first number as H, second number as Acceptor," << endl;
      cout << "  Assuming numbers in file represent indices 0..N-1" << endl;
      vector<int> hList;
      vector<int> accList;
      readPairsFile( filename, hList, accList, "index0" );

      hbManager->disableCarefulChecking();

      for ( int i = 0; i < hList.size(); i++ ) {
        bool success;
        success = false;
        //if ( !hbManager->isPairExcluded( hList[i], accList[i] ) ) {
        hbManager->addConstraint( hList[i], accList[i], success );
	//}
        if ( !success ) {
	  cout << "  Note: pair " << hList[i] << " " << accList[i] << " (index 0..N-1 format) " <<
	    "was rejected" << endl;
        }
      }
}

void TextFileInput::readHydrophobicsFile_index0( string filename, PHManager *phManager ) {
      cout << "Reading hydrophobic pairs from file " << filename << endl;
      cout << "  Interpreting first number as atom1, second number as atom2," << endl;
      cout << "  Assuming numbers in file represent indices 0..N-1" << endl;
      vector<int> p1List;
      vector<int> p2List;
      readPairsFile( filename, p1List, p2List, "index0" );

      phManager->disableCarefulChecking();

      for ( int i = 0; i < p1List.size(); i++ ) {
        bool success;
        success = false;
        //if ( !phManager->isPairExcluded( p1List[i], p2List[i] ) ) {
        phManager->addConstraint( p1List[i], p2List[i], success );
	//}
        if ( !success ) {
	  cout << "Note: pair " << p1List[i] << " " << p2List[i] << " (index 0..N-1 format) " <<
	    "was rejected" << endl;
        }
      }
}

void TextFileInput::readPairsFile(
  std::string filename,
  vector<int>& firstatom,
  vector<int>& secondatom,
  string format ) {

  bool subtract1 = false;
  if ( format == "index1" ) subtract1 = true;
  else if ( format != "index0" ) {
    cout << "TextFileInput readPairsFile error: format string must be index0 or index1" << endl;
    exit(0);
  }

  ifstream pairfile( filename.c_str(), ios::in );
  if (!pairfile) {
    cout << "Could not open file: " << filename << endl;
    exit(1);
  }
  string currentline;
  int p1;
  int p2;
  while (!pairfile.eof()) {
    getline( pairfile, currentline );
    if ( currentline.size() < 3 ) continue;
    stringstream ss;
    ss.str( currentline );
    ss >> p1;
    ss >> p2;
    if ( !ss ) {
      cout << "Error reading file " << filename << endl;
      exit(0);
    }

    ss.str("");
    ss.clear();

    if ( subtract1 ) {
      //we must convert the index from 1..N format to 0..N-1 format
      firstatom.push_back( p1-1 );
      secondatom.push_back( p2-1 );
    }
    else {
      //assume index is 0..N-1 format
      firstatom.push_back( p1 );
      secondatom.push_back( p2 );
    }
  }

}

NeighborTable *TextFileInput::buildCovalentNeighborTable_FromPrmtop( const AmberPrmtop& prmtop )
{
  int N = prmtop.natom;
  NeighborTable *covalentNeighborTable = new NeighborTable( N );
  for ( int i=0; i<N; i++ ) {
    int nNeighbors = prmtop.neighborTable[i].size();
    for ( int j=0; j<nNeighbors; j++ ) {
      covalentNeighborTable->insert( i, prmtop.neighborTable[i][j].first );
    }
  }
  covalentNeighborTable->commit();
  return covalentNeighborTable;
}

PairTypeCutoffs *TextFileInput::buildPairTypeCutoffs(
    string atomTypesAssignmentFilename,
    string pairTypeCutoffsFilename,
    int nAtoms,
    double scaleFactor ) {

  PairTypeCutoffs *pairTypeCutoffs = new PairTypeCutoffs;
  string currentline;
  stringstream ss;

  ifstream pairTypeCutoffsFile( pairTypeCutoffsFilename.c_str(), ios::in );
  if ( !pairTypeCutoffsFile.good() ) {
    cout << "Could not open " << pairTypeCutoffsFilename << endl;
    exit(0);
  }

  //read first line, which lists the number of atom types
  getline( pairTypeCutoffsFile, currentline );
  ss.str( currentline );
  int nAtomTypes;
  ss >> nAtomTypes;
  ss.str("");
  ss.clear();

  //read in the default radius of each atom type.
  //Lines should list the atomtype and the radius
  int count = 0;
  int atomtype;
  double r;
  vector<double> radii( nAtomTypes );
  while ( !pairTypeCutoffsFile.eof() && count < nAtomTypes ) {
    getline( pairTypeCutoffsFile, currentline );
    if ( currentline.size() < 3 ) continue;
    ss.str( currentline );
    ss >> atomtype;
    ss >> r;

    if ( !ss || atomtype != count ) {
      cout << "Error reading radii from " << pairTypeCutoffsFilename << endl;
      exit(0);
    }
    ss.str("");
    ss.clear();

    radii[atomtype] = r*scaleFactor;
    count++;
  }

  if ( count != nAtomTypes ) {
    cout << "Error reading radii from " << pairTypeCutoffsFilename << endl;
    exit(0);
  }
  pairTypeCutoffs->setTypeRadii( radii );

  //continue reading the rest of the lines,
  //which list pair type cutoffs that override the radius-based cutoffs.
  //Lines contain atomtype1, atomtype2, cutoff
  int atomtype1;
  int atomtype2;
  double cutoff;
  while (!pairTypeCutoffsFile.eof()) {
    getline( pairTypeCutoffsFile, currentline );
    if ( currentline.size() < 5 ) continue;
    ss.str( currentline );
    ss >> atomtype1;
    ss >> atomtype2;
    ss >> cutoff;
    ss.str("");
    ss.clear();

    pairTypeCutoffs->setPairTypeCutoff( atomtype1, atomtype2, cutoff*scaleFactor );
  }
  pairTypeCutoffsFile.close();

  vector<int> atomtypes;
  buildAtomTypes( atomTypesAssignmentFilename, atomtypes );

  pairTypeCutoffs->mapAtomsToTypes( atomtypes );
  pairTypeCutoffs->check();

  return pairTypeCutoffs;
}

void TextFileInput::buildAtomTypes( string filename, vector<int>& atomtypes ) {
  ifstream atomTypesAssignmentFile( filename.c_str(), ios::in );
  if ( !atomTypesAssignmentFile.good() ) {
    cout << "Could not open " << filename << endl;
    exit(0);
  }
  string currentline;
  stringstream ss;
  int atomIndex;
  int atomType;
  int expectedAtomIndex = 0;
  while (!atomTypesAssignmentFile.eof()) {
    getline( atomTypesAssignmentFile, currentline );
    if ( currentline.size() < 3 ) continue;
    ss.str( currentline );
    ss >> atomIndex;
    ss >> atomType;
    ss.str("");
    ss.clear();

    if ( atomIndex != expectedAtomIndex ) {
      cout << "Error reading atomtypes file: Expected atom index " << expectedAtomIndex << endl;
      exit(0);
    }
    atomtypes.push_back( atomType );
    expectedAtomIndex++;
  }
  atomTypesAssignmentFile.close();

}

TargetEnergy *TextFileInput::buildTargetEnergy(
    string filename,
    RigidUnitSystem *rigidUnitSystem ) {

  TargetEnergy *targetEnergy = new TargetEnergy( rigidUnitSystem );
  ifstream targetFile( filename.c_str(), ios::in );
  if ( !targetFile.good() ) {
    cout << "Could not open " << filename << endl;
    exit(0);
  }
  string currentline;
  stringstream ss;
  int atomIndex;
  double x;
  double y;
  double z;
  while (!targetFile.eof()) {
    getline( targetFile, currentline );
    if ( currentline.size() < 7 ) continue;
    ss.str( currentline );
    ss >> atomIndex;
    ss >> x;
    ss >> y;
    ss >> z;
    ss.str("");
    ss.clear();

    targetEnergy->addTargetPoint( atomIndex, Vec3(x,y,z) );
  }
  targetEnergy->update();
  return targetEnergy;
}

map<int,int>* TextFileInput::buildAtomMap( string mapFilename ) {
  map< int, int >* indexmap = new map<int,int>;

  //read in map file
  if ( mapFilename != "" ) {
    cout << "Reading map file " << mapFilename << endl;
    cout << "  (Interpreting first number as initial state, second number as target state)" << endl;
    ifstream mapFile( mapFilename.c_str(), ios::in );
    if ( !mapFile.good() ) {
      cout << "Error opening file " << mapFilename << endl;
      exit(0);
    }

    string currentline;
    int p1;
    int p2;
    while (!mapFile.eof()) {
      getline( mapFile, currentline );
      if ( currentline.size() < 3 ) continue;
      stringstream ss;
      ss.str( currentline );
      ss >> p1;
      ss >> p2;
      ss.str("");
      ss.clear();

      (*indexmap)[p1] = p2;
    }
    mapFile.close();
  }
  return indexmap;
}

TargetMap* TextFileInput::buildTargetMap( string mapFilename ) {
  TargetMap* targetMap = new TargetMap;

  //read in map file
  if ( mapFilename != "" ) {
    cout << "Reading map file " << mapFilename << endl;
    cout << "  (Interpreting first number as initial state, second number as target state)" << endl;
    ifstream mapFile( mapFilename.c_str(), ios::in );
    if ( !mapFile.good() ) {
      cout << "Error opening file " << mapFilename << endl;
      exit(0);
    }

    string currentline;
    int p1;
    int p2;
    while (!mapFile.eof()) {
      getline( mapFile, currentline );
      if ( currentline.size() < 3 ) continue;
      stringstream ss;
      ss.str( currentline );
      ss >> p1;
      ss >> p2;
      ss.str("");
      ss.clear();

      targetMap->insert( p1, p2 );
    }
    mapFile.close();
  }
  return targetMap;
}

void TextFileInput::readAtomSet( std::string filename, std::set<int>& atomset ) {
  //read file
  if ( filename != "" ) {
    cout << "Reading file " << filename << endl;
    ifstream infile( filename.c_str(), ios::in );
    if ( !infile.good() ) {
      cout << "Error opening file " << filename << endl;
      exit(0);
    }

    string currentline;
    int p1;
    while (!infile.eof()) {
      getline( infile, currentline );
      if ( currentline.size() < 1 ) continue;
      stringstream ss;
      ss.str( currentline );
      ss >> p1;
      if ( ss ) atomset.insert( p1 );
      ss.str("");
      ss.clear();

    }
    infile.close();
  }

}
