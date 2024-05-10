/*
 * Targeter.cpp
 *
 *  Created on: Aug 21, 2009
 *      Author: dwfarrel
 */

#include "Targeter.h"
#include "ObjectBuilder_FromSettings.h"
#include "ProteinInfo.h"
#include "NeighborTable.h"
#include "RigidUnitSystem.h"
#include "ConstraintEnforcingPotential_Targeting.h"
#include "Output.h"
#include "OutputFiles.h"
#include "Settings.h"
#include "mt19937ar.h"
#include "Timer.h"
#include "Tools.h"
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

Targeter::Targeter( const Settings& settings ) :
  dynamicConstraintRecorder( NULL )
{

  //get data on the initial structure
  cout << "Loading Initial Structure Data" << endl;
  ObjectBuilder_FromSettings* objBuilder = new ObjectBuilder_FromSettings( settings );
  prot = objBuilder->getProteinInfo();
  nt = objBuilder->getCovalentNeighborTable();
  sys = objBuilder->getRigidUnitSystem();
  delete objBuilder;

  cep = new ConstraintEnforcingPotential_Targeting( sys, *prot, *nt, settings );
  if ( settings.outputconstraintlists ) cep->writeConstraints();

  cep->minim->setToleranceCondition( "maxPreconditionedGradComponent", 0.005 );
  cep->minim->setNminimizationSteps( 200 );

  cout << "Total " << prot->natoms() << " atoms " << endl;
  cout << "Targeting " << cep->targmap->src2targ().size() << " atoms" << endl;

  //initialize targeting variables to default/initial values
  retryCount = 0;
  retryLimit = 5;
  doDynamicConstraints = settings.targDynamicConstraints;
  randomPertFreq = 0;
  iter = 0;
  state = FORWARD;
  maxIter = settings.perturbRelax.Nsteps >= 0 ?
    settings.perturbRelax.Nsteps : numeric_limits<int>::max();
  delta0 = ( settings.targDelta < numeric_limits<double>::epsilon() ) ?
    0.1 : settings.targDelta;
  acceptableViolation = 0.20;
  violationTriggerExtraMinimization = 0.10;
  doReplace = false;
  doCollectivePert = false;

  t = new TargetingTools( sys, cep, nt, prot );

  //setup random perturbation sizes
  double translationalPert_Angstroms =
    settings.perturbRelax.randomCenterPerturbationSize > numeric_limits<double>::epsilon() ?
    settings.perturbRelax.randomCenterPerturbationSize : 1.00;
  double rotationalPert_Angstroms =
    settings.perturbRelax.randomRotorPerturbationSize > numeric_limits<double>::epsilon() ?
      settings.perturbRelax.randomRotorPerturbationSize : 10.0; //5.00;
  double rotationalPert_Radians = 2.0*M_PI/3.0; //M_PI/2.0;

  cout << "Setting up Random Perturbation: " << endl;
  cout << "  translational maximum distance: " << translationalPert_Angstroms << endl;
  cout << "  maximum rotational arc-length: " << rotationalPert_Angstroms << endl;
  cout << "  maximum rotational angle (deg): " << rotationalPert_Radians*180.0/M_PI << endl;
  t->setupRandomPert( translationalPert_Angstroms, rotationalPert_Angstroms, rotationalPert_Radians );

  //perform initial fit to target
  t->setupGlobalMotionRemover_Targeting();
  t->removeGlobalMotion();

  //setup the linear reaction coordinate direction,
  //and output the initial 2-D pathway projection
  cep->targetEnergy->setDirectionOfReactionCoordinate();

  //if anything is breakable, of if user selected dynamics constraints,
  //setup the constraint remover for over-stretched bonds.
  doConstraintRemover_Cutoff =
    ( settings.commonConstraintHandling == "breakable" ||
      settings.noncommonConstraintHandling == "breakable" ||
      doDynamicConstraints );

  if ( doConstraintRemover_Cutoff ) {
    t->setupConstraintRemover_Cutoff( 0.02, 0.02, 0.01 );
  }

  //only set up the random constraint remover if the dynamic constraint setting is true
  if ( doDynamicConstraints ) {
    t->setupConstraintRemover_Random( 0.10, 0.10, 0.20 );
  }

  if ( doDynamicConstraints && cep->hb && settings.outputconstraintlists ) {
    dynamicConstraintRecorder = new DynamicConstraintRecorder<Targeter>( "dynamic_hb.txt", cep->hb, this );
  }

  //if ( doDynamicConstraints ) {
  //  t->setupForbidList( 1, this );
  //}

  if ( settings.targ == "random" ) {
    randomPertFreq = 1;
    //t->setupCollectivePert();
    //doCollectivePert = true;
  }

  if ( doReplace ) t->setupReplace( t->rmsd(), 1.0 );

  //setup backtracking
  doBacktrack = false;
  flag_noMoreBacktrack = false;
  backtrackTarget = 0;
  momentumCount = 0;
  momentumLimit = 0;
  backtrackDeltaInitial = 0;
  backtrackDelta = 0;
  backtrackDeltaScaleFactor = 0;
  backtrackDeltaLimit = 0;
  backtrackAcceptableMaxDev = 0;
  backtrackAcceptableRMSD = 0;
  savedRMSD = 0;
  if ( settings.doBacktrack ) {
    doBacktrack = true;
    flag_noMoreBacktrack = false;
    backtrackTarget = 0;
    momentumCount = 0;
    momentumLimit = 200; //settings
    backtrackDeltaInitial = 1.0; //settings
    backtrackDelta = backtrackDeltaInitial;
    backtrackDeltaScaleFactor = 2.0; //settings
    backtrackDeltaLimit = 30; //settings
    backtrackAcceptableMaxDev = 5.0; //settings
    backtrackAcceptableRMSD = 0.5; //settings
    savedRMSD = numeric_limits<double>::max();
  }

  if ( doBacktrack ) {
    t->setupMomentumPert( 0.05, 0.05 );
  }

  //setup output objects
  output = new Output( settings.output, sys );
  outputFiles = new OutputFiles( sys );
  outputFiles->setupRMSDToTargetOutput( this );
  output->registerObserver( outputFiles );

}

Targeter::~Targeter() {
  delete output;
  delete outputFiles;
  delete prot;
  delete nt;
  delete sys;
  delete cep;
  delete t;
  delete dynamicConstraintRecorder;
}

double Targeter::getRMSDtoTarget() const { return cep->targetEnergy->getRMSDtoTarget(); }

void Targeter::outputStructure() {
  output->notifyStructureReady();
}

void Targeter::run() {
  state = FORWARD;
  iter = 0;

  Timer timer;
  timer.start();

  cout << cep->generateColumnHeaderString() << endl;
  cout << "Commencing Run" << endl;

  //initialize RMSD constraint
  t->setRMSDconstraint( t->rmsd() );
  t->outputSummaryLine();

  //output initial structure
  outputStructure();

  while ( state != STOP && iter <= maxIter ) {
    iter++;
    cout << "\nCycle " << iter << endl;

    //do something depending on current state
    switch( state ) {
    case FORWARD:
      forwardStep();
      break;
    case BACK:
      backwardsStep();
      break;
    case MOMENTUM:
      momentumStep();
      break;
    case STOP:
      break;
    }

  }

  timer.stop();
  cout << static_cast<double>(timer.getutime())/1000.0 << " s" << endl;

}

void Targeter::doExtraMinimizationIfNeeded() {
  if ( cep->overlapEnergy && cep->overlapEnergy->mismatch() > violationTriggerExtraMinimization ||
       cep->sharedPointsEnergy && cep->sharedPointsEnergy->mismatch() > violationTriggerExtraMinimization ||
       cep->overridingMinDist && cep->overridingMinDist->worstDistanceViolation() > violationTriggerExtraMinimization ) {

    cout << "  Structure has large violations: Trying extra minimization cycles" << endl;
    t->extraminimize();
    t->outputSummaryLine();
  }
}

bool Targeter::structureAcceptable() {
  if ( cep->overlapEnergy && cep->overlapEnergy->mismatch() > acceptableViolation ) return false;
  if ( cep->sharedPointsEnergy && cep->sharedPointsEnergy->mismatch() > acceptableViolation ) return false;
  if ( cep->overridingMinDist && cep->overridingMinDist->worstDistanceViolation() > acceptableViolation ) return false;
  return true;
}

void Targeter::switchToForward() {
  state = FORWARD;
  retryCount = 0;
  cep->targetEnergy->setForwards();
}

void Targeter::forwardStep() {
  cout << "  Forward Step" << endl;
  if ( retryCount ) cout << "  Retry number " << retryCount << endl;

  if ( retryCount == 0 ) {
    if ( doReplace ) t->replace();
    if ( doDynamicConstraints ) {
      t->addConstraints();
    }
    if ( doDynamicConstraints || doConstraintRemover_Cutoff ) {
      t->removeConstraints_Cutoff();
    }
    if ( doDynamicConstraints ) {
      t->removeConstraints_Random();
    }
    t->saveForRevert(); //if problem pert is on, we must save for revert
    t->adjustRMSDconstraint( -delta0 );
    //cep->targetEnergy->setRMSDconstraint( 99 ); // this line is to deactivate the RMSD constraint
    //if ( t->rmsd() < cep->targetEnergy->getRMSDconstraint() ) cep->targetEnergy->setRMSDconstraint( t->rmsd() ); //this line is to ratchet the constraint forward
  }

  if ( randomPertFreq && iter%randomPertFreq == 0 || retryCount > 0 ) {
    t->randomPert();
    t->problemPert();
  }
  else {
    t->minimize();
  }

  t->removeGlobalMotion();
  t->outputSummaryLine();

  doExtraMinimizationIfNeeded();

  bool good = structureAcceptable();
  if ( good ) {
    outputStructure();
  }

  //choose next state
  if ( cep->targetEnergy->getRMSDconstraint() < numeric_limits<double>::epsilon() ) {
    cout << "RMSD constraint reached 0.0" << endl;
    switchToStop();
    return;
  }

  if ( good ) {
    retryCount = 0;
    return; //stay in current state with retry count set to 0.
  }

  //If we make it here, the structure was not acceptable.
  cout << "  Structure not acceptable." << endl;
  t->revertAll();

  if ( retryCount++ < retryLimit ) return; //stay in current state with increased retry count.

  //The only way to continue is if backtracking mode is enabled.
  //In backtracking mode, we only attempt backtracking
  //if the structure is not close enough to the target.
  if ( doBacktrack && !flag_noMoreBacktrack &&
       ( cep->targetEnergy->calcMaxDeviation() > backtrackAcceptableMaxDev ||
         t->rmsd() > backtrackAcceptableRMSD ) ) {
    switchToBacktracking();
  }
  else switchToStop();
}

void Targeter::switchToBacktracking() {
  //if we get here, we got stuck going forwards.
  //Set a backtrack target and prepare for backwards state.
  state = BACK;
  retryCount = 0;

  if ( t->rmsd() < savedRMSD - 1.0 ) {
    backtrackDelta = backtrackDeltaInitial;
    savedRMSD = t->rmsd();
  }
  else {
    backtrackDelta *= backtrackDeltaScaleFactor;
    if ( backtrackDelta > backtrackDeltaLimit ) {
      backtrackDelta = backtrackDeltaLimit;
      cout << "  Backtracking: Reached backtrack limiting RMSD, so no more backtracking after this." << endl;
      flag_noMoreBacktrack = true;
    }
  }
  backtrackTarget = t->rmsd() + backtrackDelta;
  cout << "  Backtracking: Setting backwards targeting, with backtrackTarget RMSD " << backtrackTarget << endl;
  cout << "  Letting go" << endl;
  cep->targetEnergy->disable();
  t->minimize();
  //cep->removeAllBreakableConstraints();
  t->minimize();
  t->removeGlobalMotion();
  t->setRMSDconstraint( t->rmsd() );
  t->outputSummaryLine();
  cep->targetEnergy->setBackwards();
}

void Targeter::backwardsStep() {
  cout << "  Backwards Step" << endl;
  if ( retryCount ) cout << "  Retry number " << retryCount << endl;

  if ( retryCount == 0 ) {
    t->saveForRevert();
    t->adjustRMSDconstraint( delta0 );
  }

  if ( randomPertFreq && iter%randomPertFreq == 0 || retryCount > 0 ) {
    t->randomPert();
    t->problemPert();
  }
  else {
    t->minimize();
  }

  t->removeGlobalMotion();
  t->removeConstraints_Cutoff();
  t->outputSummaryLine();

  doExtraMinimizationIfNeeded();

  bool good = structureAcceptable();
  if ( good ) {
    outputStructure();
  }

  //choose next state
  if ( good ){
    if ( t->rmsd() >= backtrackTarget ) {
      cout << "Reached backtrack target." << endl;
      switchToMomentum();
      return;
    }
    else {
      retryCount = 0;
      return; //stay in current state, set retry count to 0
    }
  }

  //If we make it here, the structure was not acceptable.
  cout << "  Structure not acceptable." << endl;
  t->revertAll();

  //if the retry number is low enough keep retrying.
  if ( retryCount++ < retryLimit ) return; //stay in current state, with retryCount increased by 1
  else {
    //Time to reverse directions, do momentum exploration for a while, and then go forwards again.
    cout << "  Backtracking: Backwards run got stuck before reaching the backtrack target" << endl;
    cout << "                There will be no more backtracking after this." << endl;
    flag_noMoreBacktrack = true;
    switchToMomentum();
  }
}

void Targeter::switchToMomentum() {
  state = MOMENTUM;
  retryCount = 0;
  momentumCount = 0;
  cout << "  Letting go" << endl;
  cep->targetEnergy->disable();
  t->minimize();
  t->removeGlobalMotion();
  t->setRMSDconstraint( t->rmsd() );
  t->outputSummaryLine();
  cep->targetEnergy->setForwards();
  t->resetMomentumPert();
  cout << "  Backtracking: Switching to Momentum. Performing " << momentumLimit << " momentum cycles with RMSD constraint." << endl;
}

void Targeter::momentumStep() {
  //do momentum cycle
  cout << "  Momentum Exploration" << endl;
  t->momentumPert();
  t->outputSummaryLine();

  doExtraMinimizationIfNeeded();

  momentumCount++;
  outputStructure();

  if ( momentumCount > momentumLimit ) {
    cout << "  Backtracking: Switching to Forward Steps" << endl;
    switchToForward();
  }
}

void Targeter::switchToStop() {
  cout << "Stopping" << endl;
  state = STOP;
}

