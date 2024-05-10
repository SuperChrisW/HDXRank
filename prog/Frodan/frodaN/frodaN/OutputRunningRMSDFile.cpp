/*
 * OutputRunningRMSDFile.cpp
 *
 *  Created on: Jan 16, 2009
 *      Author: dwfarrel
 */

#include "OutputRunningRMSDFile.h"
#include "ProteinInfo.h"
#include "RigidUnitSystem.h"
#include <iomanip>
#include <sstream>
#include <cmath>

using namespace std;

OutputRunningRMSDFile::OutputRunningRMSDFile(
  const ProteinInfo *prot,
  const RigidUnitSystem* sys_,
  string fileprefix_ ) :
    rigidUnitSystem( sys_ ),
    outputPeriod( 0 ),
    iteration( 0 ),
    sum_of_msd( 0 ),
    fileprefix( fileprefix_ ) {

  int natoms = prot->natoms();
  for ( int i = 0; i < natoms; i++ ) {
    if ( prot->atom(i).name() == "CA" || prot->atom(i).name() == "P" ) {
      caList.push_back( i );
      resiList.push_back( prot->atom(i).resi().index() + 1 );
    }
  }

  int nca = caList.size();
  caSumSquareDeviation.resize( nca, 0 );
  initialPositions.resize( nca );
  for ( int i = 0; i < nca; i++ ) {
    int ca = caList[i];
    initialPositions[i] = rigidUnitSystem->meanPositions( ca );
  }

  rmsdfile.open("RunningRMSD.txt", ios::out);
  rmsdfile << "Iteration# Current Cumulative\n";
  rmsdfile << showpoint << fixed << setprecision(3);
  rmsdfile << 0 << " " << 0 << " " << 0 << '\n';
}

OutputRunningRMSDFile::~OutputRunningRMSDFile() {
  rmsdfile.close();
}

void OutputRunningRMSDFile::recordRMSD() {
  iteration++;
  double sumSquareDeviation = 0;
  int nca = caList.size();
  for ( int i = 0; i < nca; i++ ) {
    int ca = caList[i];
    double caSquareDeviation =
      rigidUnitSystem->meanPositions( ca ).dist2( initialPositions[i] );
    sumSquareDeviation += caSquareDeviation;
    caSumSquareDeviation[i] += caSquareDeviation;
  }
  double msd = sumSquareDeviation/nca;
  sum_of_msd += msd;

  double rmsd = sqrt( msd );
  double cumulative_rmsd = sqrt( sum_of_msd/iteration );
  rmsdfile << iteration << " " << rmsd << " " << cumulative_rmsd << '\n';

  if ( outputPeriod && iteration % outputPeriod == 0 ) {
    stringstream filename;
    filename << fileprefix << "_resrmsd_" << setw(8) << setfill('0') << iteration << ".txt";
    ofstream rmsffile( filename.str().c_str(), ios::out );
    rmsffile << showpoint << fixed << setprecision(3);
    rmsffile << "Residue # Iteration " << iteration << '\n';
    for ( int i = 0; i < nca; i++ ) {
      double rmsf = sqrt( caSumSquareDeviation[i]/iteration );
      rmsffile << resiList[i] << " " << rmsf << '\n';
    }
    rmsffile.close();
    rmsdfile << flush;
  }

}

