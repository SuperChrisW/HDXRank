#include "AmberPrmtop.h"
#include "tclap/CmdLine.h"
#include "AmberTrajectoryInput.h"
#include "PDBTrajectoryInput.h"
#include "PDB.h"
#include "ProteinInfo.h"
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;
using namespace TCLAP;


int main( int argc, char ** argv ) {
  string prmtopFilename;
  string mdcrdFilename;
  string pdblistFilename;
  string outFilename;
  vector<int> atomIDs;

  cout << "\nmeasureDistanceTraveled\n" << endl;

  try {
    CmdLine cmdLine( "", ' ', "0.1" );

    //Initialize File Settings
    ValueArg<string> prmtopFilenameArg("","prmtop","Amber prmtop file",false,"","filename",cmdLine);
    ValueArg<string> mdcrdArg("","mdcrd","Amber trajectory file",false,"","filename",cmdLine);
    ValueArg<string> pdblistArg("","pdblist","File containing a list of PDB files",false,"","filename",cmdLine);
    ValueArg<string> outfileArg("o","outfile","output filename",false,"distanceTraveled.txt","filename",cmdLine);

    // parse command line
    cmdLine.parse( argc, argv );
    prmtopFilename = prmtopFilenameArg.getValue();
    mdcrdFilename = mdcrdArg.getValue();
    pdblistFilename = pdblistArg.getValue();
    outFilename = outfileArg.getValue();

  }
  catch (TCLAP::ArgException &e)  {// catch any exceptions
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
  }

  TrajectoryInput* traj = NULL;
  ProteinInfo* prot = NULL;

  if ( pdblistFilename != "" ) {
    vector<string> pdblist;
    ifstream pdblistfile( pdblistFilename.c_str(), ios::in );
    if ( pdblistfile.is_open() ) {
      string line;
      while ( !pdblistfile.eof() ) {
        getline( pdblistfile, line );
        if ( line.size() > 0 ) pdblist.push_back(line);
      }
    }
    if ( pdblist.size() == 0 ) {
      cout << "Could not obtain list of pdb files" << endl;
      exit(0);
    }
    PDB pdb( pdblist[0] );
    prot = new ProteinInfo( pdb );
    traj = new PDBTrajectoryInput( pdblist );
  }
  else if ( prmtopFilename != "" && mdcrdFilename != "" ) {
    AmberPrmtop prmtop( prmtopFilename );
    traj = new AmberTrajectoryInput( mdcrdFilename, prmtop );
    prot = new ProteinInfo( prmtop );
  }
  else {
    cout << "Error: must supply either pdblist, " << endl;
    cout << "       or prmtop and mdcrd" << endl;
    exit(0);
  }

  //here, we might optionally set up atomIDs to only contain the C-alphas.
  //For now, default is all atoms.
  int natoms = prot->natoms();
  atomIDs.resize( natoms );
  for ( int i = 0; i < natoms; i++ ) {
    atomIDs[i] = i;
  }

  int nSubset = atomIDs.size();
  vector<double> distTraveled( nSubset, 0 );

  cout << "measuring distances" << endl;
  //loop over snapshots to calculate mean structure
  vector<Vec3> coords1;
  vector<Vec3> coords2;
  if ( traj->readNextFrame( coords1 ) ) {
    while ( traj->readNextFrame( coords2 ) ) {
      for ( int i = 0; i < nSubset; i++ ) {
        int a = atomIDs[i];
        distTraveled[i] += sqrt( coords1[a].dist2( coords2[a] ) );
      }
      coords1 = coords2;
    }
  }

  traj->close();

  cout << "Writing distances to " << outFilename << '\n'<< endl;

  ofstream outfile( outFilename.c_str() , ios::out );
  for ( int i = 0; i < nSubset; i++ ) {
    outfile << atomIDs[i] << " " << distTraveled[i] << '\n';
  }
  outfile.close();

  delete traj;
  delete prot;

  return 0;
}
