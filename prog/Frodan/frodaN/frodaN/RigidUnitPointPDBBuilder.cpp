#include "RigidUnitPointPDBBuilder.h"
#include "PDB.h"
#include "RigidUnitSystem.h"
#include "NeighborTable.h"
#include <sstream>
#include <iostream>
#include <iomanip>

using namespace std;

RigidUnitPointPDBBuilder::RigidUnitPointPDBBuilder(
    const PDB *pdb, const RigidUnitSystem *sys )
{
  newpdb = new PDB;
  
  char chain = 'A';
  
  stringstream ss;
  string segid;
  ss << right << setw(4) << 0;
  segid = ss.str();
  LineInfo lineInfo;
  
  /*
  if ( sys->nRigidUnits() > 1000 ) {
    cout << "RigidUnitPointPDB only can handle 1000 rigid units" << endl;
    exit(0);
  }
  */
  
  int ru = 0;
  for ( size_t rup = 0; rup < sys->nRigidUnitPoints(); rup++ ) {
    if ( sys->getRUfromRUP( rup ) != ru ) {
      ru = sys->getRUfromRUP( rup );

      if ( chain == 'Z' ) chain = 'A';
      else chain++;
      //ss.clear();
      //ss.str("");
      //ss << right << setw(4) << ru;
      //segid = ss.str();

      lineInfo.clear();
      lineInfo.isTerLine = true;      
      newpdb->lookupLineInfo.push_back( lineInfo );
    }
    
    int p = sys->getPfromRUP( rup );
    
    // lookup PDB atom record p.
    PDBAtomLine atomLine = pdb->atomLines[p];
    ss.clear();
    ss.str("");
    ss << rup + 1;
    ss >> atomLine.serial;

    atomLine.chainID = chain;
    //atomLine.segID = segid;
    //atomLine.tempFactor = ru;

    newpdb->atomLines.push_back( atomLine );
    
    lineInfo.clear();
    lineInfo.isAtomLine = true;
    lineInfo.atomIndex = rup;
    
    newpdb->lookupLineInfo.push_back( lineInfo );
    newpdb->conectRecords.clear();
  }
}

RigidUnitPointPDBBuilder::~RigidUnitPointPDBBuilder()
{
}

/*
void RigidUnitPointPDBBuilder::addCONECTrecords( const NeighborTable& table ) {
  //to convert p neighbor table to rup neighbor table,
  //look at each p-neighbor pair in the table
  //look at the rigid units that each of the points belongs to.
  //if the two points have any rigid units in common,
  //the corresponding rups that correspond to the points
  //in the particular rigid unit are neighbors
  
  ConectRecord conectRecord;
  for ( size_t rup = 0; rup < ; rup++ ) {
    int p = sys->getPfromRUP( rup );
    //look up ru of p
    //look up bonded neighbors of p.  If any bonded neighbors do
    //belong to ru but are in a different residue
    
    size_t nNeigh = table[p].size();
    if ( nNeigh == 0 ) continue;
    if ( nNeigh > 4 ) {
      cout << "RigidUnitPointPDBBuilder Error: CONECT record can only take 4 neighbors" << endl;
      exit(0);
    }
    conectRecord.baseAtom = rup;
    conectRecord.neighbors.resize( nNeigh );
    for ( size_t j = 0; j < nNeigh; j++ ) {
      conectRecord.neighbors[j] = table[p][j];
    }
    newpdb->conectRecords.push_back( conectRecord );
  }
}
*/
