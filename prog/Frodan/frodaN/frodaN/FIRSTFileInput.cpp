#include "FIRSTFileInput.h"
#include "PDB.h"
#include "RigidUnitSystem.h"
#include "Vec3.h"
//#include "Ropes.h"
#include "NeighborTable.h"
#include "generateOverlappingRigidUnits.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <cstdlib>

using namespace std;

class UniqueIDtoArrayIndexConverter
{
public:
  UniqueIDtoArrayIndexConverter() {}
  virtual ~UniqueIDtoArrayIndexConverter() {}
  virtual int uniqueIDtoArrayIndex( std::string uniqueid ) const = 0;
};

class FIRSTIDtoArrayIndexConverter : public UniqueIDtoArrayIndexConverter
{
public:
  FIRSTIDtoArrayIndexConverter( int natoms ) : N(natoms) {}
  virtual ~FIRSTIDtoArrayIndexConverter() {}
  int uniqueIDtoArrayIndex( std::string uniqueid ) const;
private:
  int N;
};

class PDBIDtoArrayIndexConverter : public UniqueIDtoArrayIndexConverter
{
public:
  PDBIDtoArrayIndexConverter( const PDB &pdb );
  virtual ~PDBIDtoArrayIndexConverter() {}
  int uniqueIDtoArrayIndex( std::string uniqueid ) const ;
private:
  std::map<std::string, int> mapUniqueIDtoArrayIndex;
};

int FIRSTIDtoArrayIndexConverter::uniqueIDtoArrayIndex( std::string uniqueid ) const {
  std::stringstream ss;
  int tempInt;
  ss << uniqueid;
  if ( !(ss >> tempInt) ) {
    std::cout << "Error in FIRSTIDtoArrayIndexConverter: could not convert " << uniqueid << " to int" << std::endl;;
  }

  if ( tempInt < 1 || tempInt > N ) {
    std::cout << "Error in FIRSTIDtoArrayIndexConverter: FIRST ID " << uniqueid <<
    " is not in range 1 <= id <= " << N << std::endl;
    exit(0);
  }
  return tempInt-1;
}

PDBIDtoArrayIndexConverter::PDBIDtoArrayIndexConverter( const PDB &pdb )
{
  int N = pdb.atomLines.size();
  for ( int i=0; i<N; i++ ) {
    string PDBid = pdb.atomLines[i].serial;
    map<string,int>::iterator iter;
    iter = mapUniqueIDtoArrayIndex.find( PDBid );
    if ( iter != mapUniqueIDtoArrayIndex.end() ) {
      // if we find the pdb id already mapped to a particular atom index
      cout << "Error PDB_AtomIDtranslator: PDB ID \"" << PDBid << "\" is not unique." << endl;
      exit(0);
    }
    mapUniqueIDtoArrayIndex[PDBid] = i;
  }
}

int PDBIDtoArrayIndexConverter::uniqueIDtoArrayIndex( std::string uniqueid ) const {
  std::map<std::string, int>::const_iterator iter;
  iter = mapUniqueIDtoArrayIndex.find( uniqueid );
  if ( iter == mapUniqueIDtoArrayIndex.end() ) {
    std::cout << "Error UniqueIDtoArrayIndexConverter: looking up an id that does not exist" << std::endl;
    std::cout << "Bad ID: " << uniqueid << std::endl;
    exit(0);
  }
  return iter->second;
}

FIRSTFileInput::FIRSTFileInput( string pdbfilename, string FIRSTBondFileInterpretation ) :
  converter(NULL)
{
  internalpdb = new PDB( pdbfilename );
  pdb = internalpdb;
  setup( FIRSTBondFileInterpretation );
}

FIRSTFileInput::FIRSTFileInput( const PDB* pdb_, string FIRSTBondFileInterpretation ) :
  internalpdb(NULL),
  converter(NULL)
{
  pdb = pdb_;
  setup( FIRSTBondFileInterpretation );
}

FIRSTFileInput::~FIRSTFileInput()
{
  delete internalpdb;
  delete converter;
}

void FIRSTFileInput::setup( string FIRSTBondFileInterpretation ) {
  delete converter;

  int npoints = pdb->atomLines.size();
  initialPoints.resize( npoints );
  for ( int i = 0; i < npoints; i++ ) {
    initialPoints[i] = pdb->atomLines[i].position;
  }

  if ( FIRSTBondFileInterpretation == "pdb" ) {
    converter = new PDBIDtoArrayIndexConverter( *pdb );
  }
  else if ( FIRSTBondFileInterpretation == "index1" ) {
    converter = new FIRSTIDtoArrayIndexConverter( npoints );
  }
  else {
    cout << "Error: UniqueIDtoArrayIndexConverter has not been set" << endl;
    exit(0);
  }
}

/*
Ropes *FIRSTFileInput::buildRopes_FromFIRSThphobes(
    string phobesfilename,
    const RigidUnitSystem *rigidUnitSystem )
{
  if ( !pdb ) {
    cout << "Error: no Reference PDB file specified" << endl;
    exit(0);
  }

  Ropes *ropes = new Ropes( rigidUnitSystem );

  ifstream pairfile( phobesfilename.c_str(), ios::in );
  if (!pairfile) {
    cout << "Could not open file: " << phobesfilename << endl;
    exit(1);
  }
  string currentline;
  string uniqueid1;
  string uniqueid2;
  int p1;
  int p2;
  while (!pairfile.eof()) {
    getline( pairfile, currentline );
    if ( currentline.size() < 3 ) continue;
    stringstream ss;
    ss.str( currentline );
    ss >> uniqueid1;
    ss >> uniqueid2;
    ss.str("");
    ss.clear();
    p1 = converter->uniqueIDtoArrayIndex(uniqueid1);
    p2 = converter->uniqueIDtoArrayIndex(uniqueid2);
    if ( p2 < p1 ) swap( p1, p2 );
    if ( p1 == p2 ) continue;

    if ( rigidUnitSystem->doPointsBelongToSameRigidUnit( p1, p2 ) ) continue;

    double r1;
    double r2;
    double hydrophobicLength;
    if ( pdb->atomLines[p1].element == "C" ) r1 = 1.70;
    else if ( pdb->atomLines[p1].element == "S" ) r1 = 1.80;
    else {
      cout << "Error - PDB element name of hydrophobic tether atom is not C or S" << endl;
      exit(0);
    }
    if ( pdb->atomLines[p2].element == "C" ) r2 = 1.70;
    else if ( pdb->atomLines[p2].element == "S" ) r2 = 1.80;
    else {
      cout << "Error - PDB element name of hydrophobic tether atom is not C or S" << endl;
      exit(0);
    }
    hydrophobicLength = r1 + r2 + 0.5;
    double dist = sqrt( initialPoints[p1].dist2( initialPoints[p2] ) );
    if ( dist > hydrophobicLength ) {
      cout << p1 << " " << pdb->atomLines[p1].element << " " <<
              p2 << " " << pdb->atomLines[p2].element << " " << dist <<
              " > " << hydrophobicLength << endl;
    }
    ropes->insert(p1, p2, hydrophobicLength);
  }
  return ropes;
}
*/

NeighborTable* FIRSTFileInput::buildCovalentNeighborTable_FromFIRSTcov( string firstcovfilename )
{
  if ( !pdb ) {
    cout << "Error: no Reference PDB file specified" << endl;
    exit(0);
  }

  int npoints = initialPoints.size();
  NeighborTable *covalentNeighborTable = new NeighborTable( npoints );

  ifstream pairfile( firstcovfilename.c_str(), ios::in );
  if (!pairfile) {
    cout << "Could not open file: " << firstcovfilename << endl;
    exit(1);
  }
  string currentline;
  string uniqueid1;
  string uniqueid2;
  int index1;
  int index2;
  while (!pairfile.eof()) {
    getline( pairfile, currentline );
    if ( currentline.size() < 3 ) continue;
    stringstream ss;
    ss.str( currentline );
    ss >> uniqueid1;
    ss >> uniqueid2;
    ss.str("");
    ss.clear();
    index1 = converter->uniqueIDtoArrayIndex(uniqueid1);
    index2 = converter->uniqueIDtoArrayIndex(uniqueid2);
    if ( index2 < index1 ) swap( index1, index2 );
    if ( index1 == index2 ) continue;
    covalentNeighborTable->insert( index1, index2 );
  }
  covalentNeighborTable->commit();
  return covalentNeighborTable;
}

NeighborTable *FIRSTFileInput::buildHbondNeighborTable( string filename, double energyCutoff ) {

  if ( !pdb ) {
    cout << "Error: no Reference PDB file specified" << endl;
    exit(0);
  }

  int npoints = initialPoints.size();
  NeighborTable *hbondsNeighborTable = new NeighborTable( npoints );

  ifstream pairfile( filename.c_str(), ios::in );
  if (!pairfile) {
    cout << "Could not open file: " << filename << endl;
    exit(1);
  }
  string currentline;
  string uniqueid1;
  string uniqueid2;
  int index1;
  int index2;
  double bondenergy;
  while (!pairfile.eof()) {
    getline( pairfile, currentline );
    if ( currentline.size() < 3 ) continue;
    stringstream ss;
    ss.str( currentline );
    ss >> uniqueid1;
    ss >> uniqueid2;
    ss >> bondenergy;
    ss.str("");
    ss.clear();
    // In FIRST, the bond is excluded if bondenergy > energyCutoff
    // ( not '>=' ).
    // The same criteria is applied here.
    if ( bondenergy > energyCutoff ) continue;
    index1 = converter->uniqueIDtoArrayIndex(uniqueid1);
    index2 = converter->uniqueIDtoArrayIndex(uniqueid2);
    //here we ensure that index1 < index2.
    //this means that the pair is not necessarily H-A,
    //it could be A-H.
    if ( index2 < index1 ) swap( index1, index2 );
    if ( index1 == index2 ) continue;

    hbondsNeighborTable->insert( index1, index2 );
  }

  hbondsNeighborTable->commit();
  return hbondsNeighborTable;
}

vector< vector<int> > *FIRSTFileInput::buildRUtoPlist_FromFIRSTdatafile(
    std::string filename, const NeighborTable& neighborTable )
{
  vector<int> PtoRC;

  ifstream firstfile( filename.c_str(), ios::in );
  if (!firstfile) {
    cout << "Could not open file: " << filename << endl;
    exit(1);
  }

  int Natoms = neighborTable.size();

  string currentline;
  int firstnum;
  string orignum;
  int rigidclusternum;
  int nP = 0;
  while (!firstfile.eof()) {
    getline( firstfile, currentline );
    if ( currentline.size() < 60 ) continue;
    if ( currentline[0] == '#' ) continue;
    stringstream ss;
    ss.clear();
    ss.str( currentline );
    ss >> firstnum >> orignum >> rigidclusternum;
    if ( !ss ) {
      cout << "Parse error encountered in file " << filename << endl;
      exit(0);
    }
    if ( ++nP != firstnum ) {
      cout << "Error reading file " << filename
           << ": FIRST_num data expected to be 1,2,3..." << endl;
      exit(1);
    }
    if ( nP > Natoms ) {
      cout << "Error reading FIRST rigid cluster file.  Number of atoms\n";
      cout << "in file is inconsistent with the number of atoms found in\n";
      cout << "other input files." << endl;
      exit(0);
    }


    //convert rigidclusternum to C-Style index by subracting 1
    PtoRC.push_back(rigidclusternum-1);
  }
  firstfile.close();

  if ( nP != Natoms ) {
    cout << "Error reading FIRST rigid cluster file.  Number of atoms\n";
    cout << "in file is inconsistent with the number of atoms found in\n";
    cout << "other input files." << endl;
    exit(0);
  }

  vector< vector<int> > *RUtoPlist = new vector< vector<int> >;
  generateOverlappingRigidUnits( PtoRC, neighborTable, *RUtoPlist );

  return RUtoPlist;

}
