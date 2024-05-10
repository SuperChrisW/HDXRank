#include <fstream>
#include <iostream>
#include <cmath>
#include "PDB.h"
#include "tclap/CmdLine.h"
#include "PDBTrajectoryInput.h"
#include "Vec3.h"
#include "TextFileInput.h"
#include "TargetMap.h"

using namespace std;
using namespace TCLAP;

int main( int argc, char ** argv ) {
  string refpdbFilename1;
  string refpdbFilename2;
  string pdblistFilename;
  string mapFilename;

  try {
    CmdLine cmdLine( "", ' ', "0.1" );

    //Initialize File Settings
    ValueArg<string> refpdbFilename1Arg("","refpdb1","the reference PDB file 1",true,"","filename",cmdLine);
    ValueArg<string> refpdbFilename2Arg("","refpdb2","the reference PDB file 2",true,"","filename",cmdLine);
    ValueArg<string> pdblistArg("","pdblist","File containing a list of trajectory PDB files",true,"","filename",cmdLine);
    ValueArg<string> mapFilenameArg("m","map","atommap filename (index 0..N-1)",true,"","filename",cmdLine);

    // parse command line
    cmdLine.parse( argc, argv );
    refpdbFilename1 = refpdbFilename1Arg.getValue();
    refpdbFilename2 = refpdbFilename2Arg.getValue();
    pdblistFilename = pdblistArg.getValue();
    mapFilename = mapFilenameArg.getValue();
  }
  catch (TCLAP::ArgException &e)  {// catch any exceptions
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
  }

  //get the reference positions
  vector<Vec3> coords1;
  PDB* refpdb1 = new PDB( refpdbFilename1 );
  refpdb1->getPositions( coords1 );
  delete refpdb1;
  //get the reference positions
  vector<Vec3> coords2;
  PDB* refpdb2 = new PDB( refpdbFilename2 );
  refpdb2->getPositions( coords2 );
  delete refpdb2;

  //prepare list of pdbs
  vector<string> pdblist;
  if ( pdblistFilename != "" ) {
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
  }

  //get the targeted atom pairs
  TextFileInput textFileInput;
  TargetMap *targmap = textFileInput.buildTargetMap( mapFilename );

  map<int,int>::const_iterator it = targmap->src2targ().begin();
  vector<int> indices1;
  vector<int> indices2;
  for ( ; it != targmap->src2targ().end(); it++ ) {
    indices1.push_back( it->first );
    indices2.push_back( it->second );
  }
  delete targmap;

  //calc normalized difference vec
  size_t ntarg = targmap->src2targ().size();
  vector<Vec3> u( ntarg );
  if ( ntarg == 0 ) return 0;

  double norm2 = 0;
  for ( size_t i = 0; i < ntarg; i++ ) {
    int p = indices1[i];
    int q = indices2[i];
    u[i] = coords2[q] - coords1[p];
    norm2 += u[i].norm2();
  }
  double norm = sqrt(norm2);
  for ( size_t i = 0; i < ntarg; i++ ) {
    u[i] /= norm;
  }

  ofstream outfile( "proj.out", ios::out );

  //load the trajectory
  TrajectoryInput* traj = new PDBTrajectoryInput( pdblist );

  //for each frame in the trajectory, calculate projections
  while ( true ) {
    bool good;
    traj->readNextFrame( coords1, good );
    if ( !good ) break;

    double proj1 = 0;
    double msd = 0;
    double proj2 = 0;
    for ( size_t i = 0; i < ntarg; i++ ) {
      int p = indices1[i];
      int q = indices2[i];
      Vec3 pq = coords1[p] - coords2[q];
      msd += pq.norm2();
      proj1 += pq.dot( u[i] );
    }
    msd /= static_cast<double>( ntarg );
    proj1 /= sqrt( static_cast<double>( ntarg ) );
    proj2 = sqrt( msd - proj1*proj1 );

    outfile << proj1 << " " << proj2 << '\n';
  }
  delete traj;

  outfile.close();
  return 0;
}
