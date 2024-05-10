#include "PerturbRelaxCycle.h"
#include "MinimizeSystem.h"
#include "GenericMap.h"
#include "RandomRotorPerturber.h"
#include "RandomCenterPerturber.h"
#include "SymmetricPerturber.h"
#include "GlobalMotionRemover.h"
#include "MomentumPerturber.h"
#include "PhaseSpacePathLengthIntegrator.h"
#include "RMSD.h"
#include "Settings.h"
#include "RigidUnitSystem.h"
#include "ConstraintEnforcingPotential.h"
#include "ProteinInfo.h"
#include "NeighborTable.h"
#include "Output.h"
#include "EZDMap.h"
#include "PDB.h"
#include "EMStructure.h"
#include "MapPerturber.h"
#include "CorrPerturber.h"
#include "TextFileInput.h"
#include "ObjectBuilder_FromSettings.h"
#include <vector>

PerturbRelaxCycle::PerturbRelaxCycle( const Settings& settings ) :
  minim(NULL),
  map(NULL),
  randomRotorPerturber(NULL),
  randomCenterPerturber(NULL),
  symPert(NULL),
  globalMotionRemover(NULL),
  mom(NULL),
  pathLengthIntegrator(NULL),
  rmsdFromInitial(NULL),
  myEM(NULL),
  corrPert(NULL),
  mapPert(NULL)//,
{
  ObjectBuilder_FromSettings* objBuilder = new ObjectBuilder_FromSettings( settings );
  prot = objBuilder->getProteinInfo();
  nt = objBuilder->getCovalentNeighborTable();
  sys = objBuilder->getRigidUnitSystem();
  delete objBuilder;

  cep = new ConstraintEnforcingPotential( sys, *prot, *nt, settings );

  output = new Output( settings.output, sys );

  TextFileInput textFileInput;

  //Configure the Perturb/Relax Cycle.

  //initial minimization
  cep->reportConstraintViolations();
  minim = new MinimizeSystem( sys, cep );
  minim->setToleranceCondition( settings.perturbRelax.tolType, settings.perturbRelax.tol );
  minim->setNminimizationSteps( 200 );

  cout << "Performing minimization..." << flush;
  minim->minimize();
  cout << "Done." << endl;

  cep->reportConstraintViolations();

  //add minimization to the cycle
  minim->setNminimizationSteps( settings.perturbRelax.NminimizationSteps );
  addCommand_minimize( minim, &MinimizeSystem::minimize );

  if ( settings.perturbRelax.doRandomRotorPerturbation ) {
    randomRotorPerturber = new RandomRotorPerturber( sys, settings.perturbRelax.randomRotorPerturbationSize );
    addCommand_perturb( randomRotorPerturber, &RandomRotorPerturber::perturb );
  }
  if ( settings.perturbRelax.doRandomCenterPerturbation ) {
    randomCenterPerturber = new RandomCenterPerturber( sys, settings.perturbRelax.randomCenterPerturbationSize );
    addCommand_perturb( randomCenterPerturber, &RandomCenterPerturber::perturb );
  }

  if ( settings.perturbRelax.doSymmetricPerturbation) {
    symPert = new SymmetricPerturber( sys, cep->symmetryMatrices, settings.perturbRelax.symmetricPerturbationSize );
    addCommand_perturb( symPert, &SymmetricPerturber::perturb);
  }
  if ( settings.files.ezdMap != "") {
    map = new EZDMap( settings.files.ezdMap );
  }
  // if working with a correlation coefficient, set up the EMStructure object
  if (settings.perturbRelax.doCorrPerturbation) {
    PDB pdb( settings.files.pdb );
    myEM = new EMStructure(pdb,*sys,settings.perturbRelax.emResolution,
                           settings.perturbRelax.emCutoff);
  }
  if (settings.perturbRelax.doCorrPerturbation) {
    corrPert = new CorrPerturber(sys,map,myEM,settings.perturbRelax.corrPerturbationSize);
    if (settings.perturbRelax.doSymmetricPerturbation) {
      corrPert->associateMatrices(cep->symmetryMatrices);
      addCommand_perturb(corrPert, &CorrPerturber::symPerturb);
    } else {
      addCommand_perturb(corrPert, &CorrPerturber::perturb);
    }
  }
  if (settings.perturbRelax.doMapPerturbation) {
    mapPert = new MapPerturber(sys,map,settings.perturbRelax.mapPerturbationSize);
    if (settings.perturbRelax.doSymmetricPerturbation) {
      mapPert->associateMatrices(cep->symmetryMatrices);
      addCommand_perturb(mapPert, &MapPerturber::symPerturb);
    } else {
      addCommand_perturb(mapPert, &MapPerturber::perturb);
    }
  }

  if ( settings.perturbRelax.doMomentumPerturbation ) {
    mom = new MomentumPerturber( sys );
    addCommand_perturb( mom, &MomentumPerturber::perturb );
    addCommand_cycleStart( mom, &MomentumPerturber::determineDeltaQ );
    addCommand_cycleStart( mom, &MomentumPerturber::setQ1 );
  }
  if ( ( settings.perturbRelax.removeGlobalMotion ||
         settings.perturbRelax.doMomentumPerturbation )
       && !globalMotionRemover
       && settings.files.ezdMap == "" ) {
    globalMotionRemover = new GlobalMotionRemover( sys );

    if ( settings.perturbRelax.fitSelection == "CA" ) {
      //set the c-alpha mask
      vector<size_t> indexmask;
      PDB pdb( settings.files.pdb );
      size_t N = pdb.atomLines.size();
      for ( size_t i = 0; i < N; i++ ) {
        if ( pdb.atomLines[i].strippedName() == "CA" ) indexmask.push_back( i );
      }
      globalMotionRemover->setIndexMask( indexmask );
    }

    //set the current points as target, just once
    globalMotionRemover->setCurrentPointsAsTarget();

    addCommand_postMinimization( globalMotionRemover, &GlobalMotionRemover::fitCurrentPointsToTarget );
    //addCommand_cycleStart( globalMotionRemover, &GlobalMotionRemover::setCurrentPointsAsTarget );
  }

  if ( settings.perturbRelax.stopAtPhysicalTimeLimit ) {
    pathLengthIntegrator = new PhaseSpacePathLengthIntegrator(
      settings.files.pdb,
      this,
      &sys->meanPositions(),
      settings.perturbRelax.physicalTimeLimitPS );
  }

  if ( settings.perturbRelax.Nsteps > 0 ) {
    setStopAtNSteps( settings.perturbRelax.Nsteps );
  }

  rmsdFromInitial = new RMSD( sys );
  addCommand_postMinimization( this, &PerturbRelaxCycle::outputSummaryLine );
}

PerturbRelaxCycle::~PerturbRelaxCycle()
{
  delete minim;
  delete map;
  delete randomRotorPerturber;
  delete randomCenterPerturber;
  delete symPert;
  delete globalMotionRemover;
  delete mom;
  delete pathLengthIntegrator;
  delete rmsdFromInitial;
  delete myEM;
  delete corrPert;
  delete mapPert;
  delete prot;
  delete nt;
  delete cep;
  delete sys;
  delete output;
}

void PerturbRelaxCycle::outputSummaryLine() {
  cout << "Cycle " << getCycleCount() << endl;
  cout << "Minimization cycles: " << minim->getFinalStepNum() << endl;
  cout << "Violations: " << cep->generateSummaryString() << endl;
  cout << "RMSD " << rmsdFromInitial->calcRMSD() << endl;
}
