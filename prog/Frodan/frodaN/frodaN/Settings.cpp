#include "Settings.h"
#include "tclap/CmdLine.h"
#include <limits>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "ParseXmlAttributes.h"

using namespace TCLAP;
using namespace std;

class MyCommandLineOutput : public StdOutput
{
public:
  virtual void failure(CmdLineInterface& c, ArgException& e)
  {
    std::cerr << "PARSE ERROR: " << e.argId() << std::endl
              << "             " << e.error() << std::endl << std::endl;
    briefusage();
  }

  void briefusage() {
    cout << "USAGE instructions under development." << endl;
  }

  virtual void usage(CmdLineInterface& c)
  {
    briefusage();
  }
};

Settings::Settings( int argc, char **argv ) {
  string xmlfilename;

  try {
    CmdLine cmdLine( "", ' ', "none" );
    MyCommandLineOutput myCmdOut;
    cmdLine.setOutput( &myCmdOut );

    UnlabeledValueArg<string> xmlArg("xml","options filename (xml)",false,"","filename",cmdLine);

    vector<string> allowedRunValues;
    allowedRunValues.push_back( "fixedcon_old" );
    allowedRunValues.push_back( "target" );
    allowedRunValues.push_back( "dyncon" );
    allowedRunValues.push_back( "fixedcon" );
    ValuesConstraint<string> runConstraint( allowedRunValues );
    ValueArg<string> runArg("r","run","run type",false,"none",&runConstraint,cmdLine);

    vector<string> allowedTargValues;
    allowedTargValues.push_back( "direct" );
    allowedTargValues.push_back( "random" );
    ValuesConstraint<string> targConstraint( allowedTargValues );
    ValueArg<string> targArg("","targ","targ type",false,"direct",&targConstraint,cmdLine);

    vector<string> allowedConstraintHandlingValues;
    allowedConstraintHandlingValues.push_back( "fixed" );
    allowedConstraintHandlingValues.push_back( "breakable" );
    allowedConstraintHandlingValues.push_back( "remove" );
    ValuesConstraint<string> allowedConstraintHandling( allowedConstraintHandlingValues );
    ValueArg<string> commonConstraintHandlingArg("","commonConstraintHandling","handling of common constraints in targeting",false,"fixed",&allowedConstraintHandling,cmdLine);
    ValueArg<string> noncommonConstraintHandlingArg("","noncommonConstraintHandling","handling of noncommon constraints in targeting",false,"remove",&allowedConstraintHandling,cmdLine);

    //Initialize File Settings
    ValueArg<string> prmtopArg("","prmtop","Amber prmtop file",false,"","filename",cmdLine);
    ValueArg<string> restartpdbArg("","restartpdb","PDB restart coordinates file",false,"","filename",cmdLine);
    ValueArg<string> FIRSTcovArg("","cov","FIRST cov file",false,"","filename",cmdLine);
    ValueArg<string> pdbArg("","pdb","PDB file",false,"","filename",cmdLine);

    vector<string> allowedFIRSTBondFileInterpretationValues;
    allowedFIRSTBondFileInterpretationValues.push_back( "pdb" );
    allowedFIRSTBondFileInterpretationValues.push_back( "index1" );
    ValuesConstraint<string> FIRSTBondFileInterpretationConstraint( allowedFIRSTBondFileInterpretationValues );
    ValueArg<string> FIRSTBondFileInterpretationArg("","FIRSTBondFileInterpretation","Interpretation of numbers appearing in FIRST bond files.  pdb - Interpret as PDB atom IDs (can be integers, hybrid-36, or any other 5-character field as long as each atom in PDB file has a unique id), index1 (default)- Interpret as array indices 1..N.",false,"index1",&FIRSTBondFileInterpretationConstraint,cmdLine);
    ValueArg<string> FIRSThbondsArg("","hbonds","FIRST hbonds file",false,"","filename",cmdLine);
    ValueArg<string> FIRSTphobesArg("","phobes","FIRST phobes file",false,"","filename",cmdLine);
    ValueArg<string> FIRSTrcArg("","rc","FIRST rigid clusters (*_data.txt) file",false,"","filename",cmdLine);
    ValueArg<string> lessThanConstraintsFileArg("","lessThanConstraints","File listing inequality less-than constraints",false,"","filename",cmdLine);
    ValueArg<double> hbondEnergyCutoffArg("E","E","Hydrogen Bond Energy Cutoff",false,-1.0,"float",cmdLine);

    ValueArg<string> symArg("","sym","File with symmetry matrices",false,"","filename",cmdLine);
    ValueArg<string> ezdArg("","ezd","Electron density map in EZD format",false,"","filename",cmdLine);

    ValueArg<string> overrideMinDistFileArg("","overrideMinDistFile","File listing repulsive cutoffs for specific atom pairs",false,"","filename",cmdLine);
    ValueArg<string> targOverrideMinDistFileArg("","targOverrideMinDistFile","File listing repulsive cutoffs for specific atom pairs in target structure",false,"","filename",cmdLine);
    ValueArg<string> atomTypesAssignmentArg("","atomTypesAssignment","File listing atom types assignment",false,"","filename",cmdLine);
    ValueArg<string> pairTypeCutoffsArg("","pairTypeCutoffs","File listing atom types assignment",false,"","filename",cmdLine);

    ValueArg<string> atommapFileArg("","atommap","File listing atom map between initial and target structures",false,"","filename",cmdLine);
    SwitchArg targSubsetHeavyArg("","targSubsetHeavy","only apply target RMSD energy term to heavy atoms",cmdLine);
    SwitchArg doBacktrackArg("","doBacktrack","activate backtracking",cmdLine);
    ValueArg<string> targfirstcovFileArg("","targfirstcov","FIRST cov file for target structure",false,"","filename",cmdLine);
    ValueArg<string> targpdbFileArg("","targpdb","PDB file for target structure",false,"","filename",cmdLine);
    ValueArg<string> targprmtopFileArg("","targprmtop","AmberPrmtop file for target structure",false,"","filename",cmdLine);
    ValueArg<string> targAmberRestartFileArg("","targAmberRestart","Amber Restart file for target structure",false,"","filename",cmdLine);

    ValueArg<double> targDeltaArg("","targDelta","targeting - delta rmsd",false,-1,"float",cmdLine);

    //Initialize PerturbRelax Settings
    vector<string> tolType_allowedValues;
    tolType_allowedValues.push_back("maxPreconditionedGradComponent");
    ValuesConstraint<string> tolTypeConstraint( tolType_allowedValues );
    ValueArg<string> tolTypeArg("","tolType","Minimization tolerance type",false,"maxPreconditionedGradComponent",&tolTypeConstraint,cmdLine);

    vector<string> fitSelection_allowedValues;
    fitSelection_allowedValues.push_back("CA");
    fitSelection_allowedValues.push_back("all");
    ValuesConstraint<string> fitSelectionConstraint( fitSelection_allowedValues );
    ValueArg<string> fitSelectionArg("","fitSelection","atoms to use as basis for global fitting",false,"CA",&fitSelectionConstraint,cmdLine);

    ValueArg<double> tolArg("","tol","minimization tolerance",false,0.01,"float",cmdLine);
    ValueArg<int> NstepsArg("","Nsteps","number of Perturb/Relax steps",false,-1,"int",cmdLine);
    ValueArg<int> NminimizationStepsArg("","NminimizationSteps","max number of steps in each minimization",false,1000,"int",cmdLine);
    ValueArg<int> startStepArg("","startStepOffset","Begin the Perturb/Relax steps with this step number ",false,1,"int",cmdLine);
    ValueArg<double> pertCSizeArg("","pertCSize","Random Center Perturbation Size",false,0.0,"float",cmdLine);
    ValueArg<double> pertRSizeArg("","pertRSize","Random Scaled-Rotor Perturbation Size",false,0.0,"float",cmdLine);
    ValueArg<double> pertSymSizeArg("","pertSymSize","Symmetric Perturbation Size",false,0.0,"float",cmdLine);
    ValueArg<double> pertMapSizeArg("","pertMapSize","Cryo-EM map gradient perturbation size",false,0.0,"float",cmdLine);
    ValueArg<double> pertCorrSizeArg("","pertCorrSize","Cryo-EM correlation perturbation size",false,0.0,"float",cmdLine);
    ValueArg<double> emResolutionArg("","emRes","Cryo-EM resolution",false,-1.0,"float",cmdLine);
    ValueArg<double> emCutoffArg("","emCutoff","Fraction of cryo-EM density to include in distance cutoff",false,0,"float",cmdLine);
    ValueArg<double> physicalTimeLimitPSArg("","physicalTimeLimitPS","Physical Time Limit in picoseconds",false,0.0,"float",cmdLine);

    SwitchArg momentumPertArg("","momentumPert","turn on momentum perturbation",cmdLine);
    SwitchArg noRemoveGlobalMotionArg("","noRemoveGlobalMotion","do not remove global translations/rotations",cmdLine);

    SwitchArg ramaOffArg("","ramaOff","disable Rama constraints",cmdLine);
    SwitchArg torsionOffArg("","torsionOff","disable torsion constraints",cmdLine);
    SwitchArg hbondOffArg("","hbondOff","disable Hydrogen bond constraints",cmdLine);
    SwitchArg phobicOffArg("","phobicOff","disable hydrophobic constraints",cmdLine);

    //Initialize output settings
    ValueArg<int> outputConformerPeriodArg("","outputConf","How often to output a conformation",false,0,"int",cmdLine);
    SwitchArg outputRMSDFilesArg("","outputRMSDFiles","output Running RMSD file and residue RMSD files",cmdLine);
    SwitchArg amberOutputArg("","amberOutput","output AMBER trajectory file instead of PDBs",cmdLine);
    ValueArg<int> newTrajPeriodArg("","newTrajPeriod","How often to begin a new AMBER trajectory file",false,0,"int",cmdLine);
    ValueArg<double> outputConfAtRMSDFromLastArg("","outRMSDFromLast","Output conformation when a certain RMSD is reached relative to last conformation",false,0,"float",cmdLine);

    SwitchArg unbreakableInitialConArg("","unbreakableInitialCon","make initial constraints unbreakable during dynamic-constraint runs",cmdLine);
    ValueArg<int> switchToBreakArg("","switchToBreak","Iteration number (0 or greater) at which to switch to the break-constraints phase",false,-1,"int",cmdLine);
    ValueArg<int> switchToAddArg("","switchToAdd","Iteration number (0 or greater) at which to switch to the add-constraints phase",false,-1,"int",cmdLine);
    ValueArg<int> switchOffArg("","switchOff","Iteration number (0 or greater) at which to stop adding/removing constraints",false,-1,"int",cmdLine);
    ValueArg<int> switchBackToIterZeroArg("","switchBackToIterZero","Iteration number (0 or greater) at which to switch back to iteration number zero (to restart the dynamic constraints cycle)",false,-1,"int",cmdLine);

    // Initialize other settings
    ValueArg<string> seedArg("","seed","seed for Random Number Generator",false,"","hex (8 digits)",cmdLine);

    vector<string> allowedRepulsionTypeArgs;
    allowedRepulsionTypeArgs.push_back( "froda" );
    allowedRepulsionTypeArgs.push_back( "pairtype" );
    allowedRepulsionTypeArgs.push_back( "none" );
    ValuesConstraint<string> repulsionConstraint( allowedRepulsionTypeArgs );
    ValueArg<string> repulsionTypeArg("","repulsion","Type of Repulsion",false,"pairtype",&repulsionConstraint,cmdLine);
    ValueArg<double> pairCutoffScaleFactorArg("","pairCutoffScaleFactor","Pair Cutoff Scale Factor",false,1.0,"float",cmdLine);

    SwitchArg autoOverrideMinDistOffArg("","autoOverrideMinDistOff","turn off automatic detection of problem min dist constraints to override",cmdLine);
    SwitchArg redetermineConstraintsOffArg("","redetermineConstraintsOff","turn off automatic initial re-determination of constraints",cmdLine);

    ValueArg<int> targRandomPertFreqArg("","targRandomPertFreq","random perturbation frequency in targeting",false,0,"int",cmdLine);
    SwitchArg graduallyBreakNoncommonArg("","graduallyBreakNoncommon","break non-common constraints gradually",cmdLine);
    SwitchArg targDynamicConstraintsArg("","targDynamicConstraints","turn on dynamic add/removal of (breakable) constraints during targeting",cmdLine);

    SwitchArg outputconstraintlistsArg("","outputconstraintlists","turn on output of constraint lists",cmdLine);
    ValueArg<string> targatomtypesArg("","targatomtypes","targ atom types",false,"","filename",cmdLine);
    ValueArg<string> targfirstrcArg("","targfirstrc","targ first rc",false,"","filename",cmdLine);

    // parse command line
    cmdLine.parse( argc, argv );

    runtype = runArg.getValue();
    targ = targArg.getValue();

    files.targatomtypes = targatomtypesArg.getValue();
    files.targfirstrc = targfirstrcArg.getValue();
    commonConstraintHandling = commonConstraintHandlingArg.getValue();
    noncommonConstraintHandling = noncommonConstraintHandlingArg.getValue();

    doAutoOverrideMinDist = !autoOverrideMinDistOffArg.getValue();
    doRedetermineConstraints = !redetermineConstraintsOffArg.getValue();

    outputconstraintlists = outputconstraintlistsArg.getValue();

    // Retrieve File Settings
    files.hbondEnergyCutoff = hbondEnergyCutoffArg.getValue();
    files.prmtop = prmtopArg.getValue();
    files.firstcov = FIRSTcovArg.getValue();
    files.pdb = pdbArg.getValue();
    files.firsthbonds = FIRSThbondsArg.getValue();
    files.firstphobes = FIRSTphobesArg.getValue();
    files.firstrc = FIRSTrcArg.getValue();
    files.restartpdb = restartpdbArg.getValue();
    files.symMatrices = symArg.getValue();
    files.ezdMap = ezdArg.getValue();
    files.FIRSTBondFileInterpretation = FIRSTBondFileInterpretationArg.getValue();
    files.overrideMinDistFile = overrideMinDistFileArg.getValue();
    files.targOverrideMinDistFile = targOverrideMinDistFileArg.getValue();
    files.atomTypesAssignment = atomTypesAssignmentArg.getValue();
    files.pairTypeCutoffs = pairTypeCutoffsArg.getValue();
    files.atommap = atommapFileArg.getValue();
    files.targfirstcov = targfirstcovFileArg.getValue();
    files.targpdb = targpdbFileArg.getValue();
    files.targprmtop = targprmtopFileArg.getValue();
    files.targAmberRestart = targAmberRestartFileArg.getValue();

    //retrieve PerturbRelax settings
    perturbRelax.tolType = tolTypeArg.getValue();
    perturbRelax.tol = tolArg.getValue();
    perturbRelax.Nsteps = NstepsArg.getValue();
    perturbRelax.NminimizationSteps = NminimizationStepsArg.getValue();
    perturbRelax.startStep = startStepArg.getValue();
    perturbRelax.randomCenterPerturbationSize = pertCSizeArg.getValue();
    perturbRelax.doRandomCenterPerturbation =
      ( abs(perturbRelax.randomCenterPerturbationSize) > numeric_limits<double>::epsilon() );
    perturbRelax.randomRotorPerturbationSize = pertRSizeArg.getValue();
    perturbRelax.doRandomRotorPerturbation =
      ( abs(perturbRelax.randomRotorPerturbationSize) > numeric_limits<double>::epsilon() );
    perturbRelax.doMomentumPerturbation = momentumPertArg.getValue();
    perturbRelax.removeGlobalMotion = !noRemoveGlobalMotionArg.getValue();
    perturbRelax.physicalTimeLimitPS = physicalTimeLimitPSArg.getValue();
    perturbRelax.stopAtPhysicalTimeLimit = perturbRelax.physicalTimeLimitPS > numeric_limits<double>::epsilon();
    perturbRelax.symmetricPerturbationSize = pertSymSizeArg.getValue();
    perturbRelax.doSymmetricPerturbation =
      ( abs(perturbRelax.symmetricPerturbationSize) > numeric_limits<double>::epsilon() );
    perturbRelax.mapPerturbationSize = pertMapSizeArg.getValue();
    perturbRelax.doMapPerturbation =
      ( abs(perturbRelax.mapPerturbationSize) > numeric_limits<double>::epsilon() );
    perturbRelax.corrPerturbationSize = pertCorrSizeArg.getValue();
    perturbRelax.doCorrPerturbation =
      ( abs(perturbRelax.corrPerturbationSize) > numeric_limits<double>::epsilon() );
    perturbRelax.emResolution = emResolutionArg.getValue();
    perturbRelax.emCutoff = emCutoffArg.getValue();

    perturbRelax.fitSelection = fitSelectionArg.getValue();

    output.outputConformerPeriod = outputConformerPeriodArg.getValue();
    output.amberOutput = amberOutputArg.getValue();
    output.newTrajPeriod = newTrajPeriodArg.getValue();
    output.outputRMSDFiles = outputRMSDFilesArg.getValue();
    output.outputConfAtRMSDFromLast = outputConfAtRMSDFromLastArg.getValue();
    output.pdbfilename = files.pdb;
    output.prmtopfilename = files.prmtop;

    //obtain random seed
    //first try using the command line
    seedstring = seedArg.getValue();

    //retrieve other settings
    repulsionType = repulsionTypeArg.getValue();

    //
    energy.doOverlapEnergy = repulsionType != "none";
    energy.doSymmetryEnergy = files.symMatrices != "";
    energy.doRama = !ramaOffArg.getValue();
    energy.doTorsion = !torsionOffArg.getValue();
    energy.doHbond = !hbondOffArg.getValue();
    energy.doHydrophobic = !phobicOffArg.getValue();

    unbreakableInitialCon = unbreakableInitialConArg.getValue();
    switchToBreak = switchToBreakArg.getValue();
    switchToAdd = switchToAddArg.getValue();
    switchOff = switchOffArg.getValue();
    switchBackToIterZero = switchBackToIterZeroArg.getValue();
    targSubsetHeavy = targSubsetHeavyArg.getValue();
    targDelta = targDeltaArg.getValue();
    doBacktrack = doBacktrackArg.getValue();
    targDynamicConstraints = targDynamicConstraintsArg.getValue();

    pairCutoffScaleFactor = pairCutoffScaleFactorArg.getValue();

    //set default settings values (only for those settings not listed in the command line parser)
    energy.doHbondAngles = true;
    files.targsubset = "";
    files.breakableconstraints_atomlist = "";
    files.hbonds_index0 = "";
    files.hydrophobics_index0 = "";
    perturbRelax.momentumScaleFactor = 1.0;
    mergeConstraintsForTargeting = true;

    xmlfilename = xmlArg.getValue();

  }
  catch (TCLAP::ArgException &e)  {// catch any exceptions
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
  }

  try {
    if ( xmlfilename != "" ) {
      parseXml( xmlfilename );
    }
  }
  catch ( ticpp::Exception &e ) {
    std::cerr << e.what() << std::endl;
    throw e;
  }

  if ( (files.hbonds_index0 != "" || files.hydrophobics_index0 != "") && doRedetermineConstraints ) {
    cout << "Error: if using settings hydrophobics_index0 or hbonds_index0, you must also set 'redetermine' to 'off'" << endl;
    exit(0);
  }
  

  std::istringstream iss( seedstring );
  if ( seedstring=="auto" || !(iss >> hex >> seed) ) {
    //if the command line did not provide a valid seed, try reading from
    //the random "device" /dev/random (linux only)
    ifstream r( "/dev/random", ios::in | ios::binary );
    r.read( reinterpret_cast<char*>( &seed ), sizeof( unsigned long ) );

    //if that didn't work, quit.
    if ( !r ) {
      cout << "Could not obtain random seed from /dev/random," << endl;
      cout << "Must provide seed on command line with option --seed, or in xml input file" << endl;
      exit(0);
    }
    r.close();
  }

  sanityCheck();

}

Settings::~Settings()
{
}

void Settings::parseXml( string filename ) {
  cout << "Parsing settings from XML file " << filename << endl;
  ticpp::Document xmldoc( filename );
  xmldoc.LoadFile();
  ticpp::Element* root = xmldoc.FirstChildElement();
  ticpp::Iterator<ticpp::Element> e;
  cout << root->Value() << endl;
  for ( e = e.begin( root ); e != e.end(); e++ ) {
    string name;
    name = e->Value();
    cout << " xml element " << name << endl;
    if ( name == "modelfiles" ) {
      parseXMLSettings_modelfiles( e.Get() );
    }
    else if ( name == "targetfiles" ) {
      parseXMLSettings_targetfiles( e.Get() );
    }
    else if ( name == "runtarget" ) {
      parseXMLSettings_runtarget( e.Get() );
    }
    else if ( name == "runfixedconstraints" ) {
      parseXMLSettings_runfixedconstraints( e.Get() );
    }
    else if ( name == "random" ) {
      parseXMLSettings_random( e.Get() );
    }
    else if ( name == "output" ) {
      parseXMLSettings_output( e.Get() );
    }
    else if ( name == "constraints" ) {
      parseXMLSettings_constraints( e.Get() );
    }
    else {
      cout << "  (skipped)" << endl;
    }
  }

}

void Settings::parseXMLSettings_targetfiles( ticpp::Element* e ) {
  ParseXmlAttributes p;
  p.add( "pdb", &files.targpdb );
  p.add( "cov", &files.targfirstcov );
  p.add( "rigidunits", &files.targfirstrc );
  p.add( "atomtypes", &files.targatomtypes );
  p.add( "map", &files.atommap );

  p.parse( e );
}

void Settings::parseXMLSettings_modelfiles( ticpp::Element* e ) {
  ParseXmlAttributes p;
  p.add( "pdb", &files.pdb );
  p.add( "cov", &files.firstcov );
  p.add( "rigidunits", &files.firstrc );
  p.add( "atomtypes", &files.atomTypesAssignment );
  p.add( "prmtop", &files.prmtop );

  p.parse( e );
  output.pdbfilename = files.pdb;
  output.prmtopfilename = files.prmtop;
}

void Settings::parseXMLSettings_constraints( ticpp::Element* e ) {
  ParseXmlAttributes p;
  p.add( "pairtypecutoffs", &files.pairTypeCutoffs );
  p.add( "paircutoffscalefactor", &pairCutoffScaleFactor );
  p.add( "hbondenergycutoff", &files.hbondEnergyCutoff );
  p.add( "hbondangleconstraints", &energy.doHbondAngles );
  p.add( "hbond", &energy.doHbond );
  p.add( "hbondfile_index0", &files.hbonds_index0 );
  p.add( "hydrophobic", &energy.doHydrophobic );
  p.add( "hydrophobicsfile_index0", &files.hydrophobics_index0 );
  p.add( "rama", &energy.doRama );
  p.add( "torsion", &energy.doTorsion );
  p.add( "repulsion", &repulsionType );
  p.add( "override", &doAutoOverrideMinDist );
  p.add( "redetermine", &doRedetermineConstraints );

  p.parse( e );
}

void Settings::parseXMLSettings_output( ticpp::Element* e ) {
  ParseXmlAttributes p;
  p.add( "outrmsdfromlast", &output.outputConfAtRMSDFromLast );
  p.add( "outevery", &output.outputConformerPeriod );
  p.add( "amberoutput", &output.amberOutput );
  p.add( "outputrmsdfiles", &output.outputRMSDFiles );
  p.add( "outputconstraintlists", &outputconstraintlists );
  p.parse( e );
}
void Settings::parseXMLSettings_random( ticpp::Element* e ) {
  ParseXmlAttributes p;
  p.add( "seed", &seedstring );
  p.parse( e );
}

void Settings::parseXMLSettings_runtarget( ticpp::Element* e ) {
  runtype = "target";

  ParseXmlAttributes p;
  p.add( "targheavyatomsonly", &targSubsetHeavy );
  p.add( "dynamicconstraints", &targDynamicConstraints );
  p.add( "delta", &targDelta );
  p.add( "targ", &targ );
  p.add( "dobacktrack", &doBacktrack );
  p.add( "mergeconstraints", &mergeConstraintsForTargeting );
  p.add( "commonconstrainthandling", &commonConstraintHandling );
  p.add( "noncommonconstrainthandling", &noncommonConstraintHandling );
  p.add( "nsteps", &perturbRelax.Nsteps );
  p.add( "targsubset", &files.targsubset );
  p.add( "breakableconstraints_atomlist", &files.breakableconstraints_atomlist );
  p.parse( e );
}

void Settings::parseXMLSettings_runfixedconstraints( ticpp::Element* e ) {
  runtype = "fixedcon";

  ParseXmlAttributes p;
  p.add( "domomentumpert", &perturbRelax.doMomentumPerturbation );
  p.add( "momentumscalefactor", &perturbRelax.momentumScaleFactor );
  p.add( "randomperttranslation", &perturbRelax.randomCenterPerturbationSize );
  p.add( "randompertrotation", &perturbRelax.randomRotorPerturbationSize );
  p.add( "nsteps", &perturbRelax.Nsteps );
  p.parse( e );
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Checks for consistency among command-line options.  So far, I've only
//   included the options involved in cryo-EM fitting and symmetric
//   perturbations; add more checks as the need arises.
////////////////////////////////////////////////////////////////////////////////
void Settings::sanityCheck() const {
  if (perturbRelax.doMapPerturbation && perturbRelax.doCorrPerturbation) {
    cerr << "ERROR: Cryo-EM fitting should be done with either the map\n"
         << "gradient perturbation or the correlation gradient perturbation,\n"
         << "but not both.  Don't use the pertMapSize and pertCorrSize flags\n"
         << "together.\n";
    exit(1);
  }
  if (perturbRelax.doCorrPerturbation && perturbRelax.emResolution == -1.0) {
    cerr << "ERROR: Cryo-EM fitting using correlation gradient perturbation\n"
         << "requires a specified map resolution.  Always use the edRes and\n"
         << "pertCorrSize flags together.\n";
    exit(1);
  }
  if (perturbRelax.doCorrPerturbation && perturbRelax.emResolution <= 0.0) {
    cerr << "ERROR: Argument of edRes flag must be > 0.0.\n";
    exit(1);
  }
  if (perturbRelax.emResolution != -1.0 && !perturbRelax.doCorrPerturbation) {
    cerr << "ERROR: The emRes flag is only needed in when fitting cryo-EM\n"
         << "maps using the pertCorrSize flag.\n";
    exit(1);
  }
  if (perturbRelax.doCorrPerturbation && (perturbRelax.emCutoff <= 0.0 ||
                                          perturbRelax.emCutoff >= 1.0)) {
    cerr << "ERROR: Cryo-EM fitting using correlation gradient perturbation\n"
         << " requires a cutoff value between 0 and 1.  Use the emCutoff flag.\n";
    exit(1);
  }
  if (perturbRelax.emCutoff > 0.0 && !perturbRelax.doCorrPerturbation) {
    cerr << "ERROR: emCutoff is only needed when doing correlation gradient\n"
         << "fitting using pertCorrSize.\n";
    exit(1);
  }
  if (perturbRelax.doMapPerturbation && files.ezdMap == "") {
    cerr << "ERROR: Cannot fit cryo-EM map without specifying a map file;\n"
         << "use the pertMapSize and ezd flags together.\n";
    exit(1);
  }
  if (perturbRelax.doCorrPerturbation && files.ezdMap == "") {
    cerr << "ERROR: Cannot fit cryo-EM map without specifying a map file;\n"
         << "use the pertCorrSize and ezd flags together.\n";
    exit(1);
  }
  if (files.ezdMap != "" && !perturbRelax.doCorrPerturbation &&
                            !perturbRelax.doMapPerturbation) {
    cerr << "ERROR: The ezd flag must be used together with either\n"
         << "pertMapSize or pertCorrSize.\n";
    exit(1);
  }
  if (perturbRelax.doSymmetricPerturbation && files.symMatrices == "") {
    cerr << "ERROR: Cannot do symmetric perturbations without matrix input\n"
         << "file!  The sym and pertSymSize flags should always be used\n"
         << "together.\n";
    exit(1);
  }
  if (files.symMatrices != "" && !perturbRelax.doSymmetricPerturbation) {
    cerr << "ERROR: Symmetry matrices are only needed when symmetric\n"
         << "perturbation is enabled.  Always use the sym and pertSymSize\n"
         << "flags together.\n";
    exit(1);
  }
  return;
}
