/*
 * DynamicConstraints.cpp
 *
 *  Created on: Jun 11, 2009
 *      Author: dwfarrel
 */

#include "DynamicConstraints.h"
#include "RigidUnitSystem.h"
#include "ConstraintEnforcingPotential.h"
#include "ProteinInfo.h"
#include "Timer.h"
#include "Settings.h"
#include "ObjectBuilder_FromSettings.h"
#include "Output.h"
#include "OutputFiles.h"
#include "Tools.h"
#include <iostream>
#include <iomanip>
#include <algorithm>

using namespace std;

DynamicConstraints::DynamicConstraints( const Settings& settings ) :
  cycle( 0 )
{
  ObjectBuilder_FromSettings* objBuilder = new ObjectBuilder_FromSettings( settings );
  prot = objBuilder->getProteinInfo();
  nt = objBuilder->getCovalentNeighborTable();
  sys = objBuilder->getRigidUnitSystem();
  delete objBuilder;

  cep = new ConstraintEnforcingPotential( sys, *prot, *nt, settings );
  if ( settings.outputconstraintlists ) cep->writeConstraints();
  cep->minim->setToleranceCondition( "maxPreconditionedGradComponent", 0.005 );
  cep->minim->setNminimizationSteps( 200 );

  //Note.  If we're going to keep a fixed set, we can
  //load them into the dynamic
  //hbonds object as unbreakable constraints.

  bbhb = cep->bbhb;
  hb = cep->hb;
  ph = cep->ph;
  hbManager = cep->hbManager;
  phManager = cep->phManager;

  //Making unbreakable constraints here is now disabled...
  //bool breakable = !settings.unbreakableInitialCon;

  dyn = settings.runtype == "dyncon";

  if ( !dyn ) {
    delete hbManager;
    hbManager = cep->hbManager = NULL;
    delete phManager;
    phManager = cep->phManager = NULL;
  }

  t = new Tools( sys, cep, nt, prot );

  //setup global motion remover
  t->setupGlobalMotionRemover_ReferenceInitialState();
  t->removeGlobalMotion();

  doMomentumPert = settings.perturbRelax.doMomentumPerturbation;
  doRandomPert = !doMomentumPert;

  if ( doMomentumPert ) {
    cout << "Setting up Momentum Perturbation" << endl;
    cout << "  momentum scale factor: " << settings.perturbRelax.momentumScaleFactor << endl;
    
    // Determine the random perturbation sizes that will go along with the momentum perturbation.
    // If the user did not specify any values, use the defaults below, which 
    // are small (0.05), as this small size is appropriate for use with the momentum perturbation.
    double translationalPert_Angstroms =
      settings.perturbRelax.randomCenterPerturbationSize > numeric_limits<double>::epsilon() ?
      settings.perturbRelax.randomCenterPerturbationSize : 0.05;
    double rotationalPert_Angstroms =
      settings.perturbRelax.randomRotorPerturbationSize > numeric_limits<double>::epsilon() ?
        settings.perturbRelax.randomRotorPerturbationSize : 0.05;
    double rotationalPert_Radians = 2.0*M_PI/3.0;

    cout << "Setting up Random Perturbation: " << endl;
    cout << "  translational maximum distance: " << translationalPert_Angstroms << endl;
    cout << "  maximum rotational arc-length: " << rotationalPert_Angstroms << endl;
    cout << "  maximum rotational angle (deg): " << rotationalPert_Radians*180.0/M_PI << endl;

    t->setupMomentumPert( translationalPert_Angstroms, rotationalPert_Angstroms, settings.perturbRelax.momentumScaleFactor );
  }
  else if ( doRandomPert ) {
    //In the Settings object, default values of perturbation size is 0.
    //If the default value is set, we will automatically set up "big random steps".
    //If there are user-entered values in the perturbation sizes, we will use those.
    //setup random perturbation sizes
    double translationalPert_Angstroms =
      settings.perturbRelax.randomCenterPerturbationSize > numeric_limits<double>::epsilon() ?
      settings.perturbRelax.randomCenterPerturbationSize : 1.00;
    double rotationalPert_Angstroms =
      settings.perturbRelax.randomRotorPerturbationSize > numeric_limits<double>::epsilon() ?
        settings.perturbRelax.randomRotorPerturbationSize : 10.00;
    double rotationalPert_Radians = 2.0*M_PI/3.0;

    cout << "Setting up Random Perturbation: " << endl;
    cout << "  translational maximum distance: " << translationalPert_Angstroms << endl;
    cout << "  maximum rotational arc-length: " << rotationalPert_Angstroms << endl;
    cout << "  maximum rotational angle (deg): " << rotationalPert_Radians*180.0/M_PI << endl;
    t->setupRandomPert( translationalPert_Angstroms, rotationalPert_Angstroms, rotationalPert_Radians );
    //t->setupCollectivePert();
  }

  output = new Output( settings.output, sys );
  outputFiles = NULL;
  if ( settings.output.outputRMSDFiles ) {
    outputFiles = new OutputFiles( sys );
    outputFiles->setupRunningRMSD( prot, settings.output.pdbfilename );
    output->registerObserver( outputFiles );
  }

  switchToBreak = settings.switchToBreak;
  switchToAdd = settings.switchToAdd;
  switchOff = settings.switchOff;
  switchBackToIterZero = settings.switchBackToIterZero;
}

DynamicConstraints::~DynamicConstraints() {
  delete t;
  delete cep;
  delete sys;
  delete nt;
  delete prot;
  delete output;
  delete outputFiles;
}

void DynamicConstraints::run( int Ncycles ) {
  //main loop
  Timer timer;
  timer.start();

  cout << cep->generateColumnHeaderString() << endl;
  if ( dyn ) outputDynamicConstraintsColumnHeaders();

  cout << "  Commencing Run" << endl;

  bool doAdd = false;
  int i = 0;
  for ( cycle = 0; cycle < Ncycles; cycle++ ) {
    cout << "Cycle " << cycle << endl;

    if ( dyn ) {
      //do this one first
      if ( i == switchBackToIterZero ) i = 0;

      if ( i == switchToBreak ) {
        cout << "  Switching To Break Phase" << endl;
        cep->removeAllBreakableConstraints();
        doAdd = false;
      }
      if ( i == switchToAdd ) {
        cout << "  Switching To Add Phase" << endl;
        doAdd = true;
      }
      if ( i == switchOff ) {
        cout << "  Switching Off Adding/Removing of Constraints" << endl;
        doAdd = false;
      }
      //when we reach here, the state has been set

      if ( doAdd ) addConstraints();
    }

    perturbSystemAndRestoreContraints();
    t->outputSummaryLine();
    if ( dyn ) outputDynamicConstraintsColumns();
    output->notifyStructureReady();
    i++;
  }
  timer.stop();
  cout << static_cast<double>(timer.getutime())/1000.0 << " s" << endl;
}

void DynamicConstraints::perturbSystemAndRestoreContraints() {
  if ( doMomentumPert ) {
    t->momentumPert();
  }
  if ( doRandomPert ) {
    t->saveForRevert();
    //t->collectivePert();
    t->randomPert();
    t->problemPert();
    t->removeGlobalMotion();
  }
}

void DynamicConstraints::addConstraints() {
  if ( hbManager ) {
    hbManager->tighten();
    hbManager->findnew();
  }
  if ( phManager ) {
    phManager->tighten();
    phManager->findnew();
  }
  cep->notifyTermsChanged();
}

void DynamicConstraints::outputDynamicConstraintsColumnHeaders() {
  int i = 1;
  if ( dyn )
    cout << "Dynamic Constraints Columns:\n";
  if ( dyn && hbManager && hb && bbhb ) {
    cout << i++ << " Number of Hbonds Added (before cycle)\n";
    cout << i++ << " Number of Hbonds Removed (before cycle)\n";
    cout << i++ << " Number of Hbonds Total (before cycle)\n";
  }
  if ( dyn && phManager && ph ) {
    cout << i++ << " Number of Hydrophobics Added (before cycle)\n";
    cout << i++ << " Number of Hydrophobics Removed (before cycle)\n";
    cout << i++ << " Number of Hydrophobics Total (before cycle)\n";
  }
  cout << endl;
}

void DynamicConstraints::outputDynamicConstraintsColumns() {
  cout << fixed << setprecision(5) << setfill(' ');

  if ( dyn && hbManager && hb && bbhb ) {
    cout << " " << setw(8) << hbManager->nAdded();
    cout << " " << setw(8) << hbManager->nRemoved();
    cout << " " << setw(8) << bbhb->size() + hb->size();
    hbManager->resetCounts();
  }
  if ( dyn && phManager && ph ) {
    cout << " " << setw(8) << phManager->nAdded();
    cout << " " << setw(8) << phManager->nRemoved();
    cout << " " << setw(8) << ph->size();
    phManager->resetCounts();
  }
  cout << endl;

}
