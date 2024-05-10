/*
#include "OutputRigidUnits.h"
#include "PerturbRelaxCycle.h"
#include "PDB.h"
#include "MinimizeSystem.h"
#include "RigidUnitSystem.h"
#include "RigidUnitPointPDBBuilder.h"
#include <sstream>
#include <iostream>
#include <iomanip>

using namespace std;

OutputRigidUnits::OutputRigidUnits(
    string inputpdbfilename,
    const RigidUnitSystem *rigidUnitSystem_,
    const PerturbRelaxCycle *cycle_ ) :
  pdb(),
  rigidUnitSystem(rigidUnitSystem_),
  cycle( cycle_ ),
  countConf(0),
  minCycleCount(0)
{
  PDB inputpdb( inputpdbfilename );
  RigidUnitPointPDBBuilder rupPDBBuilder( &inputpdb, rigidUnitSystem );
  pdb = rupPDBBuilder.getPDB();
}

OutputRigidUnits::~OutputRigidUnits()
{
  delete pdb;
}

void OutputRigidUnits::write( string filename ) {
  pdb->setPositions( rigidUnitSystem->absolutePositions() );
  pdb->write( filename );
}

void OutputRigidUnits::beforePerturb() {
  stringstream ss;
  ss << "iter" << setfill('0') << setw(8) << ++countConf << ".A" << ".pdb";
  write(ss.str());
  minCycleCount = 0;
}

void OutputRigidUnits::afterPerturb() {
  stringstream ss;
  ss << "iter" << setfill('0') << setw(8) << countConf << ".B" << ".pdb";
  write(ss.str());
}

void OutputRigidUnits::afterMinCycle() {
  stringstream ss;
  ss << "iter" << setfill('0') << setw(8) << countConf << ".C." <<
        setw(5) << ++minCycleCount << ".pdb";
  write(ss.str());
}

void OutputRigidUnits::finishedConformer() {
  stringstream ss;
  ss << "iter" << setfill('0') << setw(8) << countConf << ".D" << ".pdb";
  write(ss.str());
}

void OutputRigidUnits::receiveNotification( Observable *obs ) {
  if ( obs == cycle ) {
    if ( cycle->getState() == 0 ) beforePerturb();
    else if ( cycle->getState() == 2 ) afterPerturb();
    else if ( cycle->getState() == 4 ) finishedConformer();
    return;
  }
  if ( obs == cycle->minim ) {
    afterMinCycle();
    return;
  }
}
*/
