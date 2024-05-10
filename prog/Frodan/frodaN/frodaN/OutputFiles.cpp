#include "OutputFiles.h"
#include "RigidUnitSystem.h"
#include "Targeter.h"
#include "OutputRunningRMSDFile.h"
#include "ProteinInfo.h"
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

OutputFiles::OutputFiles( const RigidUnitSystem *rigidUnitSystem_ ) :
  rigidUnitSystem(rigidUnitSystem_),
  snapshot( 1 ),
  targ( NULL ),
  outputRunningRMSDFile( NULL )
{
}

OutputFiles::~OutputFiles()
{
  delete outputRunningRMSDFile;
}

void OutputFiles::receiveNotification( Observable* obs ) {
  writeFiles();
  snapshot++;
}

void OutputFiles::setupRMSDToTargetOutput( const Targeter *targ_ ) {
  targ = targ_;
  filenameRMSDToTarget = "rmsdToTarget.txt";
  fileRMSDToTarget.open( filenameRMSDToTarget.c_str(), ios::out );
  fileRMSDToTarget.close();
}

void OutputFiles::setupRunningRMSD( const ProteinInfo *prot, string pdbfilename ) {
  string fileprefix = pdbfilename.substr( 0, pdbfilename.find_last_of('.') );
  outputRunningRMSDFile = new
    OutputRunningRMSDFile( prot, rigidUnitSystem, fileprefix );
  outputRunningRMSDFile->setOutputPeriod( 1 );
}

void OutputFiles::writeRMSDToTarget() {
  fileRMSDToTarget.open( filenameRMSDToTarget.c_str(), ios::app );
  fileRMSDToTarget << snapshot << " " << targ->getRMSDtoTarget() << '\n';
  fileRMSDToTarget.close();
}

void OutputFiles::writeFiles() {
  if ( targ ) {
    writeRMSDToTarget();
  }
  if ( outputRunningRMSDFile ) {
    outputRunningRMSDFile->recordRMSD();
  }
}
