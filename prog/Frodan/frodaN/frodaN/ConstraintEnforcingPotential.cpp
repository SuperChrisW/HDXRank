#include "ConstraintEnforcingPotential.h"
#include "SharedPoints.h"
#include "OverlapEnergy.h"
#include "SymmetryEnforcer.h"
#include "NeighborTable.h"
#include "HBConstraint.h"
#include "HBManager.h"
#include "PHConstraint.h"
#include "PHManager.h"
#include "SideChainTorsion.h"
#include "Rama.h"
#include "SymmetryMatrices.h"
#include "Repulsion.h"
#include "Energy.h"
#include "Gradient.h"
#include "Mismatch.h"
#include "RigidUnitSystem.h"
#include "Vec3.h"
#include "Settings.h"
#include "PairTypeCutoffs.h"
#include "FrodaCutoffsBuilder.h"
#include "TextFileInput.h"
#include "MinimizeSystem.h"
#include "DistConstraint.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cstdlib>
#include <sys/stat.h>

using namespace std;

ConstraintEnforcingPotential::ConstraintEnforcingPotential(
    RigidUnitSystem* rigidUnitSystem,
    const ProteinInfo& prot,
    const NeighborTable& nt,
    const Settings& settings ) :
      coords( rigidUnitSystem->meanPositions() ),
      sharedPointsEnergy(NULL),
      overlapEnergy(NULL),
      overridingMinDist(NULL),
      overridingSC(NULL),
      overridingRama(NULL),
      repulsion(NULL),
      pairTypeCutoffs(NULL),
      bbhb(NULL),
      hb(NULL),
      hbManager(NULL),
      ph(NULL),
      phManager(NULL),
      sideChainTorsion(NULL),
      rama(NULL),
      symmetryEnforcerEnergy(NULL),
      symmetryMatrices(NULL)
{

  //Set up Energy, Gradient, Mismatch
  energy_ = new Energy();
  gradient_ = new Gradient( rigidUnitSystem );
  mismatch_ = new Mismatch();
  rigidUnitSystem->registerObserver( energy_ );
  rigidUnitSystem->registerObserver( gradient_ );
  rigidUnitSystem->registerObserver( mismatch_ );

  //set up Shared Points Energy
  sharedPointsEnergy = new SharedPoints( rigidUnitSystem );
  addEnergyTerm( sharedPointsEnergy );
  addGradientTerm( sharedPointsEnergy );
  addMismatchTerm( sharedPointsEnergy );
  cout << "Shared Points energy term established" << endl;

  //set up Symmetry energy
  if ( settings.energy.doSymmetryEnergy && settings.files.symMatrices != "" ) {
    symmetryMatrices = new SymmetryMatrices(settings.files.symMatrices);
    symmetryEnforcerEnergy = new SymmetryEnforcer(rigidUnitSystem,
        symmetryMatrices);
    rigidUnitSystem->registerObserver( symmetryEnforcerEnergy );
    addEnergyTerm( symmetryEnforcerEnergy );
    addGradientTerm( symmetryEnforcerEnergy );
    addMismatchTerm( symmetryEnforcerEnergy );
    cout << "Symmetry energy term established" << endl;
  }

  //set up pair-type cutoffs
  TextFileInput textFileInput;
  if ( settings.files.atomTypesAssignment != ""  &&
       settings.files.pairTypeCutoffs != "" ) {

    struct stat st;
    int statcode;
    string filename = settings.files.pairTypeCutoffs;
    statcode = stat( filename.c_str(), &st );

    //if the filename is not found using the given path,
    //try using the FRODANHOME environment variable,
    //looking in the $FRODANHOME/dat directory
    if ( statcode != 0 ) {
      char* gemshome = getenv( "FRODANHOME" );
      if ( gemshome ) {
        filename = gemshome;
        char sep = '/';
        //make sure the FRODANHOME path ends in the file separation character
        if ( filename.size() > 0 && filename[filename.size()-1] != sep ) filename += sep;
        //now create the path
        filename += "dat";
        filename += sep;
        filename += settings.files.pairTypeCutoffs;
      }
    }

    pairTypeCutoffs = textFileInput.buildPairTypeCutoffs(
        settings.files.atomTypesAssignment,
        filename,
        rigidUnitSystem->nPoints(),
        settings.pairCutoffScaleFactor );
    cout << "Using pair cutoffs from file " << settings.files.pairTypeCutoffs <<
            " with scale factor " << settings.pairCutoffScaleFactor << endl;
  }
  else if ( settings.repulsionType == "froda" ) {
    FrodaCutoffsBuilder f( prot, nt, 0.0, settings.pairCutoffScaleFactor );
    pairTypeCutoffs = f.getPairTypeCutoffs();
    cout << "Using original FRODA radii with scale factor " << settings.pairCutoffScaleFactor << '\n';
    cout << "  [scale factor in original FRODA was 0.90 (notice contact) and 0.85 (cutoff)]" << endl;
  }

  //setup overlap energy, with verlet list and repulsion object
  if ( pairTypeCutoffs ) {
    repulsion = new Repulsion( rigidUnitSystem, &nt, pairTypeCutoffs );

    //setup overlap energy object
    overlapEnergy = new OverlapEnergy( rigidUnitSystem, repulsion );
    addEnergyTerm( overlapEnergy );
    addGradientTerm( overlapEnergy );
    addMismatchTerm( overlapEnergy );

    cout << "Overlap energy term established" << endl;
  }

  int natoms = prot.natoms();

  //set side chain torsion constraints
  if ( settings.energy.doTorsion ) {
    sideChainTorsion = new SideChainTorsionContainer( natoms );
    SideChainTorsionInitializer* torsionInit = new SideChainTorsionInitializer(
      prot, nt, rigidUnitSystem->meanPositions(), *sideChainTorsion );
    torsionInit->setupConstraints();
    torsionInit->excludeRigidPairs( rigidUnitSystem );
    delete torsionInit;
    if ( repulsion ) {
      for ( SideChainTorsionContainer::const_iterator it = sideChainTorsion->begin();
            it != sideChainTorsion->end(); it++ ) {
        repulsion->exclude( it->getp1(), it->getp2() ); //exclude
      }
    }
    addEnergyTerm( sideChainTorsion );
    addGradientTerm( sideChainTorsion );
    cout << "Side Chain Torsion energy term established, " << sideChainTorsion->size() << " pairs." << endl;
  }

  //setup Ramachandran constraints
  if ( settings.energy.doRama ) {
    rama = new RamaContainer( rigidUnitSystem, prot, nt );
    if ( repulsion ) {
      for ( RamaContainer::const_iterator it = rama->begin();
            it != rama->end(); it++ ) {
        repulsion->exclude( it->getp1(), it->getp2() ); //exclude
      }
    }
    addEnergyTerm( rama );
    addGradientTerm( rama );
    cout << "Ramachandran energy term established, " << rama->size() << " pairs." << endl;
  }

  overridingMinDist = new MinDistConstraintContainer( natoms );
  addEnergyTerm( overridingMinDist );
  addGradientTerm( overridingMinDist );
  overridingSC = new MinDistConstraintContainer( natoms );
  addEnergyTerm( overridingSC );
  addGradientTerm( overridingSC );
  overridingRama = new MinDistConstraintContainer( natoms );
  addEnergyTerm( overridingRama );
  addGradientTerm( overridingRama );


  //setup overriding minimum-distance constraints
  if ( settings.files.overrideMinDistFile != "" ) {
    cout << "Overriding minimum-distance constraints from file " <<
      settings.files.overrideMinDistFile << endl;
    setupOverridingMinDist( settings.files.overrideMinDistFile );
  }

  //setup hydrogen bond constraints
  if ( settings.energy.doHbond ) {
    bbhb = new BBHBContainer( natoms );
    hb = new HBContainer( natoms );
    hbManager = new HBManager( prot, nt, rigidUnitSystem->meanPositions(), *bbhb, *hb );
    if ( settings.energy.doHbondAngles )
      hbManager->enableAnglesInNewConstraints();
    else hbManager->disableAnglesInNewConstraints();
    hbManager->setEcutoff( settings.files.hbondEnergyCutoff );
    TextFileInput textFileInput;
    if ( settings.files.hbonds_index0 != "" ) {
      textFileInput.readHBondFile_index0( settings.files.hbonds_index0, hbManager );
    } 
    else hbManager->findnew();
    addEnergyTerm( bbhb );
    addEnergyTerm( hb );
    addGradientTerm( bbhb );
    addGradientTerm( hb );
    cout << "Hydrogen Bond energy term established, cutoff energy " <<
         settings.files.hbondEnergyCutoff << ", " << bbhb->size() + hb->size() << " pairs." << endl;
  }

  //setup hydrophobic constraints
  if ( settings.energy.doHydrophobic ) {
    ph = new PHContainer( natoms );
    phManager = new PHManager( prot, nt, rigidUnitSystem->meanPositions(), *ph );
    TextFileInput textFileInput;
    if ( settings.files.hydrophobics_index0 != "" ) {
      textFileInput.readHydrophobicsFile_index0( settings.files.hydrophobics_index0, phManager );
    }
    else  phManager->findnew();
    addEnergyTerm( ph );
    addGradientTerm( ph );
    cout << "Hydrophobic energy term established, " << ph->size() << " pairs." << endl;
  }

  //before any minimization, override min dist constraints to be
  //consistent with initial input coordinates.
  if ( settings.doAutoOverrideMinDist ) {
    cout << "Adding overrides for problem min dist constraints" << endl;
    overrideProblemMinDistConstraints();
  }

  //setup minimizer
  minim = new MinimizeSystem( rigidUnitSystem, this );
  minim->setToleranceCondition( "maxPreconditionedGradComponent", 0.001 );
  minim->setNminimizationSteps( 200 );

  //do initial minimization
  if ( settings.doRedetermineConstraints ) {

    cout << "Performing minimization..." << endl;
    minimize();

    int count = 0;

    //if there are some large hbond or phobic violations, adjust constraints and re-minimize
    while ( count < 3 &&
            (  bbhb && bbhb->worstDistanceViolation() > 0.02 ||
               hb && hb->worstDistanceViolation() > 0.02 ||
               ph && ph->worstDistanceViolation() > 0.02 ) ) {
      if ( hbManager && bbhb && hb ) {
        cout << "Redetermining HB constraint values" << endl;
        bbhb->clear();
        hb->clear();
        hbManager->findnew();
        cout << "Hydrogen Bond pairs " << bbhb->size() + hb->size() << endl;
      }
      if ( phManager && ph ) {
        cout << "Redetermining PH constraint values" << endl;
        ph->clear();
        phManager->findnew();
        cout << "Hydrophobic energy term, " << ph->size() << " pairs." << endl;
      }
      notifyTermsChanged();
      minimize();
      count++;
    }
  }

}

void ConstraintEnforcingPotential::writeConstraints() {
  if ( bbhb && hb ) {
    cout << "Writing h-bond constraints to files bbhb.txt and hb.txt" << endl;
    bbhb->write( "bbhb.txt" );
    hb->write( "hb.txt" );
  }
  if ( ph ) {
    cout << "Writing hydrophobic constraints to file ph.txt" << endl;
    ph->write( "ph.txt" );
  }
  if ( sideChainTorsion ) {
    cout << "Writing side chain torsion constraints to file sc.txt" << endl;
    sideChainTorsion->write( "sc.txt" );
  }
  if ( overridingMinDist && overridingMinDist->size() > 0 ) {
    cout << "Writing overriding min dist constraints to file overrideMinDist.txt" << endl;
    overridingMinDist->write("overrideMinDist.txt");
  }
  if ( overridingSC && overridingSC->size() > 0 ) {
    cout << "Appending overriding side chain constraints to file overrideMinDist.txt" << endl;
    overridingSC->append("overrideMinDist.txt");
  }
  if ( overridingRama && overridingRama->size() > 0 ) {
    cout << "Appending overriding rama constraints to file overrideMinDist.txt" << endl;
    overridingRama->append("overrideMinDist.txt");
  }
}

ConstraintEnforcingPotential::~ConstraintEnforcingPotential()
{
  delete sharedPointsEnergy;
  delete overlapEnergy;
  delete overridingMinDist;
  delete overridingSC;
  delete overridingRama;
  delete symmetryEnforcerEnergy;
  delete repulsion;
  delete symmetryMatrices;
  delete bbhb;
  delete hb;
  delete hbManager;
  delete ph;
  delete phManager;
  delete sideChainTorsion;
  delete rama;
  delete energy_;
  delete gradient_;
  delete mismatch_;
  delete pairTypeCutoffs;
}

void ConstraintEnforcingPotential::overrideConstraint( int p1, int p2, double newcutoff ) {
  //This function only DECREASES existing constraint distances.  It does not increase them.

  double k = overlapEnergy ? overlapEnergy->getk() : 10.0;

  //in order to override the constraint, we need to find it first.
  //It could be a regular repulsion constraint, or a side chain torsion constraint,
  //or a rama constraint.

  if ( overridingMinDist ) {
    MinDistConstraintContainer::iterator it;
    it = overridingMinDist->find( p1, p2 );
    if ( it != overridingMinDist->end() ) {
      if ( newcutoff < it->getCutoff() ) {
        it->setCutoff( newcutoff );
      }
      return;
    }
  }

  if ( overlapEnergy && repulsion ) {
    double cutoff;
    bool isPairExcluded;
    repulsion->getCutoff( p1, p2, isPairExcluded, cutoff );
    if ( !isPairExcluded ) {
      if ( newcutoff < cutoff ) {
        k = overlapEnergy->getk();
        repulsion->exclude( p1, p2 );
        overridingMinDist->insert( p1, p2, MinDistConstraint( &coords, k, p1, p2, newcutoff ) );
      }
      return;
    }
  }

  if ( overridingSC ) {
    MinDistConstraintContainer::iterator it;
    it = overridingSC->find( p1, p2 );
    if ( it != overridingSC->end() ) {
      if ( newcutoff < it->getCutoff() ) {
        it->setCutoff( newcutoff );
      }
      return;
    }
  }

  if ( sideChainTorsion ) {
    SideChainTorsionContainer::iterator sc = sideChainTorsion->find( p1, p2 );
    if ( sc != sideChainTorsion->end() ) {
      if ( newcutoff < sc->getCutoff() ) {
        sc->setCutoff( newcutoff );
        overridingSC->insert( p1, p2, *sc );
        sideChainTorsion->erase( sc );
      }
      return;
    }
  }

  if ( overridingRama ) {
    MinDistConstraintContainer::iterator it;
    it = overridingRama->find( p1, p2 );
    if ( it != overridingRama->end() ) {
      if ( newcutoff < it->getCutoff() ) {
        it->setCutoff( newcutoff );
      }
      return;
    }
  }

  if ( rama ) {
    RamaContainer::iterator ra = rama->find( p1, p2 );
    if ( ra != rama->end() ) {
      if ( newcutoff < ra->getCutoff() ) {
        ra->setCutoff( newcutoff );
        overridingRama->insert( p1, p2, *ra );
        rama->erase( ra );
      }
      return;
    }
  }

  //if we get here, then we did not find the constraint anywhere. ]
  //That means the constraint distance is currently zero.
  //If we override the constraint, the constraint distance will INCREASE
  //which is not allowed.
  //overridingMinDist->insert( p1, p2, MinDistConstraint( &coords, k, p1, p2, newcutoff ) );
}

void ConstraintEnforcingPotential::setupOverridingMinDist( string filename ) {
  vector<int> firstatom;
  vector<int> secondatom;
  vector<double> dist;

  TextFileInput textFileInput;
  textFileInput.readConstraintsFile(
    filename, firstatom, secondatom, dist );

  const int nconstraints = firstatom.size();
  for ( int i = 0; i < nconstraints; i++ ) {
    int p1 = firstatom[i];
    int p2 = secondatom[i];
    double newcutoff = dist[i];

    overrideConstraint( p1, p2, newcutoff );
  }

  if ( overlapEnergy ) overlapEnergy->refresh();
  notifyTermsChanged();
}

void ConstraintEnforcingPotential::overrideProblemMinDistConstraints() {
  double eps = numeric_limits<double>::epsilon();

  if ( overridingMinDist ) {
    for ( MinDistConstraintContainer::iterator it = overridingMinDist->begin();
          it != overridingMinDist->end(); it++ ) {
      double overlapDist = it->getCutoff() - it->calcDist();
      if ( overlapDist > eps ) {
        double newCutoff = it->calcDist();
        it->setCutoff( newCutoff );
      }
    }
  }

  //transfer from repulsion/overlapEnergy objects to override object
  if ( overlapEnergy && repulsion && overridingMinDist ) {
    for ( vector<MinDistConstraint>::iterator it = overlapEnergy->constraints.begin();
          it != overlapEnergy->constraints.end(); it++ ) {
      double overlapDist = it->getCutoff() - it->calcDist();
      if ( overlapDist > eps ) {
        double newCutoff = it->calcDist();
        it->setCutoff( newCutoff );
        overridingMinDist->insert( it->getp1(), it->getp2(), *it );
        repulsion->exclude( it->getp1(), it->getp2() );
      }
    }
    overlapEnergy->refresh();
  }

  //transfer from side chain object to override object, keeping side chain k constant
  if ( sideChainTorsion && overridingSC ) {
    SideChainTorsionContainer::iterator it = sideChainTorsion->begin();
    while ( it != sideChainTorsion->end() ) {
      double overlapDist = it->getCutoff() - it->calcDist();
      if ( overlapDist > eps ) {
        double newCutoff = it->calcDist();
        it->setCutoff( newCutoff );
        overridingSC->insert( it->getp1(), it->getp2(), *it );
        sideChainTorsion->erase( it );
      }
      else it++;
    }
  }

  //transfer from rama object to override object, keeping rama k constant
  if ( rama && overridingRama ) {
    RamaContainer::iterator it = rama->begin();
    while ( it != rama->end() ) {
      double overlapDist = it->getCutoff() - it->calcDist();
      if ( overlapDist > eps ) {
        double newCutoff = it->calcDist();
        it->setCutoff( newCutoff );
        overridingRama->insert( it->getp1(), it->getp2(), *it );
        rama->erase( it );
      }
      else it++;
    }
  }

  notifyTermsChanged();
}

void ConstraintEnforcingPotential::reportConstraintViolations() {
  cout << "Listing structure constraint violations (if any): " << endl;
  cout << "  (listed by internal array indices 0..N-1) " << endl;
  mismatch_->setVerbose( true );
  mismatch_->calc();
  mismatch_->setVerbose( false );
}

void ConstraintEnforcingPotential::removeAllBreakableConstraints() {
  if ( hbManager ) hbManager->breakAllBreakable();
  if ( phManager ) phManager->breakAllBreakable();
  notifyTermsChanged();
}

string ConstraintEnforcingPotential::generateColumnHeaderString() {
  ostringstream o;
  int i = 1;
  o << "Columns:\n";
  if ( sharedPointsEnergy ) o << i++ << " Worst Shared-Point Distance\n";
  if ( overlapEnergy || overridingMinDist )  o << i++ << " Worst Overlap Distance\n";
  if ( bbhb && hb ) o << i++ << " Worst H-bond Distance Violation\n";
  if ( bbhb && hb ) o << i++ << " Worst H-bond Angle Violation (degrees)\n";
  if ( ph ) o << i++ << " Worst Hydrophobic Distance Violation\n";
  if ( rama || overridingRama ) o << i++ << " Worst Rama distance violation\n";
  if ( sideChainTorsion || overridingSC ) o << i++ << " Worst side-chain torsion distance violation\n";
  o << endl;
  return o.str();
}

string ConstraintEnforcingPotential::generateSummaryString() {
  ostringstream o;
  o << fixed << setprecision(5) << setfill(' ');
  if ( sharedPointsEnergy ) o << " " << setw(8) << sharedPointsEnergy->mismatch();
  if ( overlapEnergy || overridingMinDist ) o << " " << setw(8) <<
    max( overlapEnergy ? overlapEnergy->mismatch() : 0,
      overridingMinDist ? overridingMinDist->worstDistanceViolation() : 0 );

  if ( bbhb && hb ) o << " " << setw(8) << max( bbhb->worstDistanceViolation(), hb->worstDistanceViolation() );
  if ( bbhb && hb ) o << " " << setw(8) << max( bbhb->worstAngleViolationDeg(), hb->worstAngleViolationDeg() );
  if ( ph ) o << " " << setw(8) << ph->worstDistanceViolation();
  if ( rama || overridingRama ) o << " " << setw(8) <<
    max( rama ? rama->worstDistanceViolation() : 0,
      overridingRama ? overridingRama->worstDistanceViolation() : 0 );
  if ( sideChainTorsion ) o << " " << setw(8) <<
    max( sideChainTorsion ? sideChainTorsion->worstDistanceViolation() : 0,
      overridingSC ? overridingSC->worstDistanceViolation() : 0 );
  o << endl;
  return o.str();
}

