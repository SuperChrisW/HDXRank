/* Outputs a list of numeric atom types to a file, like this:
 *   0 21
 *   1 5
 *   2 11
 *   3 21
 *   4 1
 *   ...
 *
 * The first column is the atom index (starting at 0), and the
 * second column is the atom type
 *
*/

#include "tclap/CmdLine.h"
#include "AmberPrmtop.h"
#include "NeighborTable.h"
#include "TextFileInput.h"
#include "AtomCategories.h"
#include "ProteinInfo.h"
#include "PDB.h"
#include "FIRSTFileInput.h"
#include <iostream>
#include <fstream>

using namespace TCLAP;
using namespace std;

int main( int argc, char ** argv ) {
  string pdbFilename;
  string covFilename;
  string prmtopFilename;
  string outFilename;
  string FIRSTBondFileInterpretation;
  try {
    CmdLine cmdLine( "", ' ', "0.1" );

    //Initialize File Settings
    ValueArg<string> pdbFilenameArg("","pdb","pdb filename",false,"","filename",cmdLine);
    ValueArg<string> covFilenameArg("","cov","FIRST cov bonds filename",false,"","filename",cmdLine);
    vector<string> allowedFIRSTBondFileInterpretationValues;
    allowedFIRSTBondFileInterpretationValues.push_back( "pdb" );
    allowedFIRSTBondFileInterpretationValues.push_back( "index1" );
    ValuesConstraint<string> FIRSTBondFileInterpretationConstraint( allowedFIRSTBondFileInterpretationValues );
    ValueArg<string> FIRSTBondFileInterpretationArg("","FIRSTBondFileInterpretation","Interpretation of numbers appearing in FIRST bond files.  Default: pdb - Interpret as PDB atom IDs (can be integers, hybrid-36, or any other 5-character field as long as each atom in PDB file has a unique id), index1 - Interpret as array indices 1..N.",false,"pdb",&FIRSTBondFileInterpretationConstraint,cmdLine);
    ValueArg<string> prmtopFilenameArg("","prmtop","prmtop filename",false,"","filename",cmdLine);
    ValueArg<string> outFilenameArg("o","outfile","output filename",true,"","filename",cmdLine);
    // parse command line
    cmdLine.parse( argc, argv );

    pdbFilename = pdbFilenameArg.getValue();
    covFilename = covFilenameArg.getValue();
    prmtopFilename = prmtopFilenameArg.getValue();
    outFilename = outFilenameArg.getValue();
    FIRSTBondFileInterpretation = FIRSTBondFileInterpretationArg.getValue();
  }
  catch (TCLAP::ArgException &e)  {// catch any exceptions
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
  }

  NeighborTable *neighborTable = NULL;
  PDB *pdb = NULL;
  AmberPrmtop *prmtop = NULL;
  ProteinInfo* prot = NULL;

  if ( pdbFilename != "" ) {
    cout << " Assigning atom types for " << pdbFilename << endl;
    pdb = new PDB( pdbFilename );
    prot = new ProteinInfo( *pdb );
    FIRSTFileInput firstFileInput( pdbFilename, FIRSTBondFileInterpretation );
    neighborTable = firstFileInput.buildCovalentNeighborTable_FromFIRSTcov( covFilename );
  }
  else if ( prmtopFilename != "" ) {
    cout << " Assigning atom types for " << prmtopFilename << endl;
    prmtop = new AmberPrmtop( prmtopFilename );
    prot = new ProteinInfo( *prmtop );
    TextFileInput textFileInput;
    neighborTable = textFileInput.buildCovalentNeighborTable_FromPrmtop( *prmtop );
  }
  else {
    cout << "Error: must supply either " << endl;
    cout << "  a PDB file and FIRST cov file, or" << endl;
    cout << "  an Amber7 prmtop file" << endl;
    exit(0);
  }

  AtomCategories atomCategories( neighborTable, prot );

  ofstream outfile;
  outfile.open( outFilename.c_str(), ios::out );
  int natoms = neighborTable->size();
  for ( int i = 0; i < natoms; i++ ) {
    outfile << i << " " << atomCategories.getAtomTypeForAtom( i ) << '\n';
  }
  outfile.close();

  delete neighborTable;
  delete prot;
  delete prmtop;
  delete pdb;

  return 0;
}
