#include "AmberRestartFileInput.h"
#include "AmberPrmtop.h"
#include "AmberTrajectory.h"
#include "Vec3.h"
#include "tclap/CmdLine.h"
#include "Fit.h"
#include "Rotator.h"
#include "TextFileInput.h"
#include "NeighborTable.h"
#include "AtomLookup.h"
#include "PDB.h"
#include "Swap.h"
#include "FIRSTFileInput.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <fstream>
#include <cstdlib>

using namespace std;
using namespace TCLAP;


int main( int argc, char ** argv ) {
  string prmtop1Filename;
  string prmtop2Filename;
  string restart1Filename;
  string restart2Filename;
  string pdb1Filename;
  string pdb2Filename;
  string cov1Filename;
  string FIRSTFileFormat;
  string mapFilename;
  string out1Filename;

  cout << "\nmatchState1ToState2\n" << endl;

  try {
    CmdLine cmdLine( "", ' ', "0.1" );
    vector<string> allowedFIRSTFileFormatValues;
    allowedFIRSTFileFormatValues.push_back( "pdb" );
    allowedFIRSTFileFormatValues.push_back( "index1" );
    ValuesConstraint<string> FIRSTFileFormatConstraint( allowedFIRSTFileFormatValues );
    //Initialize File Settings
    ValueArg<string> prmtop1FilenameArg("","prmtop1","Amber prmtop file for struct 1",false,"","filename",cmdLine);
    ValueArg<string> prmtop2FilenameArg("","prmtop2","Amber prmtop file for struct 2",false,"","filename",cmdLine);
    ValueArg<string> restart1FilenameArg("","restart1","Amber restart file (struct 1)",false,"","filename",cmdLine);
    ValueArg<string> restart2FilenameArg("","restart2","Amber restart file (struct 2)",false,"","filename",cmdLine);
    ValueArg<string> pdb1FilenameArg("","pdb1","PDB file (struct 1)",false,"","filename",cmdLine);
    ValueArg<string> pdb2FilenameArg("","pdb2","PDB file (struct 2)",false,"","filename",cmdLine);
    ValueArg<string> cov1FilenameArg("","FIRSTcov1","FIRST covalent bonds file (struct 1)",false,"","filename",cmdLine);
    ValueArg<string> FIRSTFileFormatArg("","FIRSTFileFormat","Interpretation of numbers appearing in FIRST bond files.  Default: pdb - Interpret as PDB atom IDs (can be integers, hybrid-36, or any other 5-character field as long as each atom in PDB file has a unique id), index1 - Interpret as array indices 1..N.",false,"pdb",&FIRSTFileFormatConstraint,cmdLine);
    ValueArg<string> mapFilenameArg("m","map","atommap filename (index 0..N-1)",false,"","filename",cmdLine);
    ValueArg<string> out1FilenameArg("o","out1","output filename (PDB or Amber restart) for modified struct 1",true,"","filename",cmdLine);

    // parse command line
    cmdLine.parse( argc, argv );
    prmtop1Filename = prmtop1FilenameArg.getValue();
    prmtop2Filename = prmtop2FilenameArg.getValue();
    restart1Filename = restart1FilenameArg.getValue();
    restart2Filename = restart2FilenameArg.getValue();
    pdb1Filename = pdb1FilenameArg.getValue();
    pdb2Filename = pdb2FilenameArg.getValue();
    cov1Filename = cov1FilenameArg.getValue();
    FIRSTFileFormat = FIRSTFileFormatArg.getValue();
    mapFilename = mapFilenameArg.getValue();
    out1Filename = out1FilenameArg.getValue();
  }
  catch (TCLAP::ArgException &e)  {// catch any exceptions
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
  }

  PDB pdb1;
  vector<Vec3> state1coords;
  NeighborTable *nt = NULL;
  vector<string> elem;
  vector<string> resnames;
  vector<int> res_begin;
  vector<string> names;
  int natoms;
  bool prmtop1_ifbox = false;

  if ( pdb1Filename != "" && cov1Filename != "" ) {
    //extract coords
    cout << "Obtaining structure 1 coordinates from PDB file " << pdb1Filename << endl;
    pdb1.read( pdb1Filename );
    pdb1.getPositions( state1coords );

    //get neighbor table
    FIRSTFileInput firstFileInput( pdb1Filename, FIRSTFileFormat );
    cout << "Obtaining structure 1 covalent bond topology from FIRST cov file " << cov1Filename << endl;
    nt = firstFileInput.buildCovalentNeighborTable_FromFIRSTcov( cov1Filename );

    //now extract other necessary data from the pdb file
    natoms = state1coords.size();
    elem.resize( natoms );
    names.resize( natoms );
    int last_resID = -1;
    char last_chainID = 0;
    char last_insertionCode = 0;
    for ( int i = 0; i < natoms; i++ ) {
      elem[i] = pdb1.atomLines[i].element;
      if ( elem[i].size() == 0 || elem[i][0] < 'A' || elem[i][0] > 'Z' ) {
        cout << "Error, element column in pdb file is missing or has bad data." << endl;
        exit(0);
      }
      names[i] = pdb1.atomLines[i].strippedName();
      int this_resID = pdb1.atomLines[i].resSeq;
      char this_chainID = pdb1.atomLines[i].chainID;
      char this_insertionCode = pdb1.atomLines[i].iCode;
      if ( this_resID != last_resID ||
           this_chainID != last_chainID ||
           this_insertionCode != last_insertionCode ) {
        last_resID = this_resID;
        last_chainID = this_chainID;
        last_insertionCode = this_insertionCode;
        res_begin.push_back( i );
        resnames.push_back( pdb1.atomLines[i].resName );
      }
    }
  }
  else if ( prmtop1Filename != "" && restart1Filename != "" ) {
    cout << "Obtaining structure 1 coordinates from Amber restart file " << restart1Filename << '\n';
    cout << "  with Amber prmtop file (version 7+) " << prmtop1Filename << endl;

    //prmtop
    AmberPrmtop prmtop( prmtop1Filename );
    prmtop1_ifbox = prmtop.ifbox;
    natoms = prmtop.natom;

    //input coordinates
    bool good;
    AmberRestartFileInput restrt;
    restrt.read( restart1Filename, prmtop, state1coords, good );
    if ( !good ) {
      cout << "Error reading restart file" << endl;
      exit(0);
    }

    TextFileInput textFileInput;
    nt = textFileInput.buildCovalentNeighborTable_FromPrmtop( prmtop );
    elem.resize( prmtop.natom );
    for ( int i = 0; i < prmtop.natom; i++ ) {
      elem[i] = prmtop.amberAtomType[i][0];
    }

    resnames = prmtop.residueLabels;
    res_begin = prmtop.residueFirstAtomIndex;
    names = prmtop.atomName;

  }
  else {
    cout << "Error: For structure 1, must supply either " << endl;
    cout << "  a PDB file and FIRST covalent bond file, or" << endl;
    cout << "  an Amber7 prmtop file and restart file" << endl;
    exit(0);
  }

  //read structure 2 coordinates
  vector<Vec3> state2coords;
  if ( pdb2Filename != "" ) {
    //extract coords
    cout << "Obtaining structure 2 coordinates from PDB file " << pdb2Filename << endl;
    PDB pdb( pdb2Filename );
    pdb.getPositions( state2coords );
  }
  else if ( prmtop2Filename != "" && restart2Filename != "" ) {
    cout << "Obtaining structure 2 coordinates from Amber restart file " << restart2Filename << '\n';
    cout << "  with Amber prmtop file (version 7+) " << prmtop2Filename << endl;

    //prmtop
    AmberPrmtop prmtop( prmtop2Filename );

    //input coordinates
    bool good;
    AmberRestartFileInput restrt;
    restrt.read( restart2Filename, prmtop, state2coords, good );
    if ( !good ) {
      cout << "Error reading restart file" << endl;
      exit(0);
    }
  }
  else {
    cout << "Error: For structure 2, must supply either " << endl;
    cout << "  a PDB file, or" << endl;
    cout << "  an Amber7 prmtop file and restart file" << endl;
    exit(0);
  }

  //read in map file
  map< int, int > indexmap;
  if ( mapFilename != "" ) {
    cout << "Reading map file " << mapFilename << endl;
    cout << "  (Interpreting first number as struct 1, second number as struct 2)" << endl;
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

      indexmap[p1] = p2;
    }
    mapFile.close();
  }

  AtomLookup alookup( *nt, elem );

  //initial fit
  vector<int> indices1;
  vector<int> indices2;
  if ( indexmap.size() == 0 ) {
    for ( int i = 0; i < natoms; i++ ) {
      indices1.push_back( i );
      indices2.push_back( i );
    }
  }
  else {
    for ( int i = 0; i < natoms; i++ ) {
      map<int,int>::const_iterator it = indexmap.find( i );
      if ( it == indexmap.end() ) continue;
      indices1.push_back( i );
      indices2.push_back( it->second );
    }
  }

  if ( indices1.size() <= 3 ) {
    cout << "Warning, about to fit structure 1 to structure 2, but\n";
    cout << "number of atoms in the atom-map is <= 3.  Fit may be\n";
    cout << "unstable." << endl;
  }

  Fit fit;
  Vec3 trans;
  Vec3 rotor;
  Vec3 centerOfRotation;
  Rotator rotator;
  fit.setSourceAbsolutePoints( state1coords, indices1 );
  fit.setTargetAbsolutePoints( state2coords, indices2 );
  fit.simpleFit();
  fit.getFitStep1_translation( trans );
  fit.getFitStep2_rotor_centerOfRotation( rotor, centerOfRotation );
  rotator.setRotor( rotor );
  for ( int i = 0; i < natoms; i++ ) {
    state1coords[i] += trans;
    rotator.rotateAboutCenter( state1coords[i], centerOfRotation, state1coords[i] );
  }

  CoordinateSwapper *coordinateSwapper = NULL;
  if ( indexmap.size() == 0 ) {
    coordinateSwapper = new CoordinateSwapper( state1coords, state2coords );
  }
  else {
    coordinateSwapper = new CoordinateSwapper( state1coords, state2coords, &indexmap );
  }
  vector<int> fittingSubset;
  SwapSet_2Fold s2;
  SwapPoint_3Fold s3;
  vector<int> H3;
  int nres = res_begin.size();
  for ( int ires = 0; ires < nres; ires++ ) {
    //cout << ires << " " << resnames[ires] << endl;
    alookup.clearfail();
    alookup.clearatoms();
    int res_end = ( ires < nres - 1 ) ? res_begin[ires+1] : natoms;
    for ( int i = res_begin[ires]; i < res_end; i++ ) {
      alookup.add( i, names[i] );
    }
    fittingSubset.clear();
    s2.clear();

    if ( resnames[ires] == "PHE" ||
         resnames[ires] == "TYR" ) {

      fittingSubset.push_back( alookup("CA") );
      fittingSubset.push_back( alookup("CB") );
      fittingSubset.push_back( alookup("CG") );

      s2.addPair( alookup("CD1"), alookup("CD2"));
      s2.addPair( alookup("CE1"), alookup("CE2"));
      s2.addPair( alookup.neighH1("CD1"), alookup.neighH1("CD2"));
      s2.addPair( alookup.neighH1("CE1"), alookup.neighH1("CE2"));

      if ( !alookup.fail() )
        coordinateSwapper->swap( fittingSubset, s2 );
    }
    else if ( resnames[ires] == "ASP" ) {

      fittingSubset.push_back( alookup("CA") );
      fittingSubset.push_back( alookup("CB") );
      fittingSubset.push_back( alookup("CG") );

      s2.addPair( alookup("OD1"), alookup("OD2"));

      if ( !alookup.fail() )
        coordinateSwapper->swap( fittingSubset, s2 );
    }
    else if ( resnames[ires] == "GLU" ) {

      fittingSubset.push_back( alookup("CB") );
      fittingSubset.push_back( alookup("CG") );
      fittingSubset.push_back( alookup("CD") );

      s2.addPair( alookup("OE1"), alookup("OE2"));

      if ( !alookup.fail() )
        coordinateSwapper->swap( fittingSubset, s2 );
    }
    else if ( resnames[ires] == "LEU" ) {

      fittingSubset.push_back( alookup("CB") );
      fittingSubset.push_back( alookup("CG") );
      fittingSubset.push_back( alookup("CD1") );

      alookup.neighH3( "CD1", H3 );
      s3.setTriple( H3 );
      if ( !alookup.fail() )
        coordinateSwapper->swap( fittingSubset, s3 );

      alookup.clearfail();
      fittingSubset.clear();
      fittingSubset.push_back( alookup("CB") );
      fittingSubset.push_back( alookup("CG") );
      fittingSubset.push_back( alookup("CD2") );

      alookup.neighH3( "CD2", H3 );
      s3.setTriple( H3 );
      if ( !alookup.fail() )
        coordinateSwapper->swap( fittingSubset, s3 );
    }
    else if ( resnames[ires] == "ILE" ) {

      fittingSubset.push_back( alookup("CB") );
      fittingSubset.push_back( alookup("CG1") );
      fittingSubset.push_back( alookup("CD1") );
      alookup.neighH3( "CD1", H3 );
      s3.setTriple( H3 );
      if ( !alookup.fail() )
        coordinateSwapper->swap( fittingSubset, s3 );

      alookup.clearfail();
      fittingSubset.clear();
      fittingSubset.push_back( alookup("CA") );
      fittingSubset.push_back( alookup("CB") );
      fittingSubset.push_back( alookup("CG2") );

      alookup.neighH3( "CG2", H3 );
      s3.setTriple( H3 );
      if ( !alookup.fail() )
        coordinateSwapper->swap( fittingSubset, s3 );
    }
    else if ( resnames[ires] == "VAL" ) {

      fittingSubset.push_back( alookup("CA") );
      fittingSubset.push_back( alookup("CB") );
      fittingSubset.push_back( alookup("CG1") );

      alookup.neighH3( "CG1", H3 );
      s3.setTriple( H3 );
      if ( !alookup.fail() )
        coordinateSwapper->swap( fittingSubset, s3 );

      alookup.clearfail();
      fittingSubset.clear();
      fittingSubset.push_back( alookup("CA") );
      fittingSubset.push_back( alookup("CB") );
      fittingSubset.push_back( alookup("CG2") );
      alookup.neighH3( "CG2", H3 );
      s3.setTriple( H3 );
      if ( !alookup.fail() )
        coordinateSwapper->swap( fittingSubset, s3 );
    }
    else if ( resnames[ires] == "ALA" ) {

      fittingSubset.push_back( alookup("C") );
      fittingSubset.push_back( alookup("CA") );
      fittingSubset.push_back( alookup("CB") );

      alookup.neighH3( "CB", H3 );
      s3.setTriple( H3 );
      if ( !alookup.fail() )
        coordinateSwapper->swap( fittingSubset, s3 );
    }
    else if ( resnames[ires] == "THR" ) {

      fittingSubset.push_back( alookup("CA") );
      fittingSubset.push_back( alookup("CB") );
      fittingSubset.push_back( alookup("CG2") );

      alookup.neighH3( "CG2", H3 );
      s3.setTriple( H3 );
      if ( !alookup.fail() )
        coordinateSwapper->swap( fittingSubset, s3 );
    }
    else if ( resnames[ires] == "LYS" ) {
      fittingSubset.push_back( alookup("CD") );
      fittingSubset.push_back( alookup("CE") );
      fittingSubset.push_back( alookup("NZ") );

      alookup.neighH3( "NZ", H3 );
      s3.setTriple( H3 );
      if ( !alookup.fail() )
        coordinateSwapper->swap( fittingSubset, s3 );
    }
    else if ( resnames[ires] == "MET" ) {

      fittingSubset.push_back( alookup("CG") );
      fittingSubset.push_back( alookup("SD") );
      fittingSubset.push_back( alookup("CE") );

      alookup.neighH3( "CE", H3 );
      s3.setTriple( H3 );
      if ( !alookup.fail() )
        coordinateSwapper->swap( fittingSubset, s3 );
    }

  }

  //global fit again, now that coordinates have been swapped.
  fit.setSourceAbsolutePoints( state1coords, indices1 );
  fit.setTargetAbsolutePoints( state2coords, indices2 );
  fit.simpleFit();
  fit.getFitStep1_translation( trans );
  fit.getFitStep2_rotor_centerOfRotation( rotor, centerOfRotation );
  rotator.setRotor( rotor );
  for ( int i = 0; i < natoms; i++ ) {
    state1coords[i] += trans;
    rotator.rotateAboutCenter( state1coords[i], centerOfRotation, state1coords[i] );
  }

  if ( pdb1Filename != "" ) {
    pdb1.setPositions( state1coords );
    pdb1.write( out1Filename );
  }
  else if ( restart1Filename != "" ) {
    AmberTrajectory amberTrajectory;
    amberTrajectory.writeRestartFile( out1Filename, prmtop1_ifbox, state1coords );
  }

  delete coordinateSwapper;

  vector<double> dist( indices1.size() );

  cout << "measuring distances" << endl;
  for ( int i = 0; i < indices1.size(); i++ ) {
    int a = indices1[i];
    int b = indices2[i];
    dist[i] += sqrt( state1coords[a].dist2( state2coords[b] ) );
  }

  string outFilename = "distances_State1ToState2.txt";
  cout << "Writing distances to " << outFilename << '\n'<< endl;

  ofstream outfile( outFilename.c_str() , ios::out );
  for ( int i = 0; i < indices1.size(); i++ ) {
    outfile << indices1[i] << " " << dist[i] << '\n';
  }
  outfile.close();

  cout << "matchState1ToState2 finished.\n" << endl;

  return 0;
}
