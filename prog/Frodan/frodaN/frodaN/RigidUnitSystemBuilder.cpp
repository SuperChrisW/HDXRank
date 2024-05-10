#include "RigidUnitSystemBuilder.h"
#include "RigidUnitSystem.h"
#include "PDB.h"
#include "AmberPrmtop.h"
#include "NeighborTable.h"
#include "Vec3.h"
#include "RigidUnitFitter.h"
#include "TextFileInput.h"
#include "FIRSTFileInput.h"
#include <iostream> 
#include <vector>
using namespace std;

RigidUnitSystem *RigidUnitSystemBuilder::build() {
  NeighborTable *covalentNeighborTable = NULL;
  NeighborTable *combinedNeighborTable = NULL;
  vector< vector<int> > *RUtoPlist = NULL;

  if ( pdbfilename == "" ) {
    cout << "Error: no PDB file specified" << endl;
    exit(0);
  }
  
  TextFileInput *textFileInput = new TextFileInput;  
  FIRSTFileInput *firstFileInput = new FIRSTFileInput( 
      pdbfilename, FIRSTBondFileInterpretation );

  int npoints = firstFileInput->getPDBAtomPositions().size();
  
  //covalentNeighborTable
  if ( prmtopfilename != "" ) {
    AmberPrmtop prmtop( prmtopfilename );
    covalentNeighborTable = textFileInput->buildCovalentNeighborTable_FromPrmtop( prmtop );
  }  
  else if ( firstcovfilename != "" ) {
    covalentNeighborTable = firstFileInput->buildCovalentNeighborTable_FromFIRSTcov( firstcovfilename );
  }  
  else {
    cout << "Error: covalent neighbor table not specified" << endl;
    exit(0);
  }  

  //combinedNeighborTable
  combinedNeighborTable = new NeighborTable( npoints );
  combinedNeighborTable->insert( *covalentNeighborTable );

  //Hbonds
  if ( firsthbondsfilename != "" ) {
    NeighborTable *hbondsNeighborTable = firstFileInput->buildHbondNeighborTable( firsthbondsfilename, hbondEnergyCutoff );
    combinedNeighborTable->insert( *hbondsNeighborTable ); 
    delete hbondsNeighborTable;
  }

  //Overlapping Rigid Units
  if ( firstrcfilename != "" ) {
    RUtoPlist = firstFileInput->buildRUtoPlist_FromFIRSTdatafile( 
        firstrcfilename, *combinedNeighborTable );
  }
  else {
    cout << "Error: FIRST rigid clusters not specified" << endl;
    exit(0);
  }

  RigidUnitSystem *rigidUnitSystem = new RigidUnitSystem( 
      *RUtoPlist,
      firstFileInput->getPDBAtomPositions() );

  if ( restartpdbfilename != "" ) {
    vector<Vec3> restartPositions;
    PDB restartpdb( restartpdbfilename );
    restartpdb.getPositions( restartPositions );
    RigidUnitFitter fitter( rigidUnitSystem );
    fitter.calcFitToMeanPoints( restartPositions );
    size_t nRU = rigidUnitSystem->nRigidUnits();
    Vec3 rotor;
    Vec3 deltaCenter;
    for ( size_t ru = 0; ru < nRU; ru++ ) {
      rotor = fitter.getFitRotation( ru );
      deltaCenter = fitter.getFitTranslation( ru );
      rigidUnitSystem->setRotor( ru, rotor );
      rigidUnitSystem->collapseRotor( ru );
      rigidUnitSystem->addToCenter( ru, deltaCenter );
    }
    rigidUnitSystem->update();
  }

  delete combinedNeighborTable;
  delete covalentNeighborTable;
  delete RUtoPlist;

  delete textFileInput;
  delete firstFileInput;

  return rigidUnitSystem;


}
