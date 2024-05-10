#include "Output.h"
#include "RigidUnitSystem.h"
#include "PDB.h"
#include "AmberTrajectory.h"
#include "AmberPrmtop.h"
#include "RMSD.h"
#include <limits>
#include <sstream>
#include <iostream>
#include <iomanip>

using namespace std;

Output::Output(
    const OutputSettings& settings,
    RigidUnitSystem *rigidUnitSystem_ ) :
  rigidUnitSystem(rigidUnitSystem_),
  pdb( NULL ),
  traj( NULL ),
  iter( 0 ),
  snapshot( 1 ),
  outputFrequency( 0 ),
  doRMSDFromLast( false ),
  rmsdTrigger( 0 ),
  rmsdFromLast( NULL )
{
  if ( settings.outputConformerPeriod > 0 ) {
    outputFrequency = settings.outputConformerPeriod;
  }
  else if ( settings.outputConfAtRMSDFromLast > numeric_limits<double>::epsilon() ) {
    doRMSDFromLast = true;
    rmsdTrigger = settings.outputConfAtRMSDFromLast;
    rmsdFromLast = new RMSD( rigidUnitSystem );
  }

  if ( settings.amberOutput && settings.prmtopfilename != "" ) {
    setupAmberTrajOutput( settings.prmtopfilename );
  }
  else if ( settings.pdbfilename != "" ){
    setupPDBOutput( settings.pdbfilename );
  }
}

Output::~Output()
{
}

void Output::notifyStructureReady() {
  bool doOutput = false;
  if ( iter == 0 ) {
    doOutput = true;
  }
  else if ( outputFrequency ) {
    doOutput = iter%outputFrequency == 0;
  }
  else if ( doRMSDFromLast ) {
    double rmsd = rmsdFromLast->calcRMSD();
    if ( rmsd > rmsdTrigger ) {
      rmsdFromLast->setCurrentPointsAsReference();
      doOutput = true;
    }
  }

  if ( doOutput ) {
    writeStructures();
    notifyObservers();
    snapshot++;
  }

  iter++;
}

void Output::writeStructures() {
  if ( pdb ) {
    writePDB();
  }
  if ( traj ) {
    writeAmberTraj();
  }
}

void Output::setupPDBOutput( string inputpdbfilename ) {
  pdb = new PDB( inputpdbfilename );
  pdbPrefix = inputpdbfilename.substr( 0, inputpdbfilename.find_last_of('.') );
}

void Output::writePDB()
{
  pdb->setPositions( rigidUnitSystem->meanPositions() );
  stringstream ss;
  ss << pdbPrefix << "_snapshot_" << setfill('0') << setw(8) << snapshot << ".pdb";
  pdb->write( ss.str() );
  cout << "  Wrote snapshot " << snapshot << endl;
}

void Output::setupAmberTrajOutput( string prmtopfilename ) {
  AmberPrmtop prmtop( prmtopfilename );
  string filename = "traj.mdcrd";
  traj = new AmberTrajectory;
  traj->initializeOutputFile( filename.c_str(), prmtop.ifbox );
}

void Output::writeAmberTraj() {
  traj->append( rigidUnitSystem->meanPositions() );
  cout << "  Wrote snapshot " << snapshot << endl;
}


