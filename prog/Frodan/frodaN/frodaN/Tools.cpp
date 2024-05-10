/*
 * Tools.cpp
 *
 *  Created on: Mar 24, 2010
 *      Author: dwfarrel
 */

#include "Tools.h"
#include "RigidUnitSystem.h"
#include "RandomRotorPerturber.h"
#include "RandomCenterPerturber.h"
#include "MomentumPerturber.h"
#include "ProteinInfo.h"
#include "NeighborTable.h"
#include "ConstraintEnforcingPotential_Targeting.h"
#include "ConstraintRemover.h"
#include "RigidUnitReplace.h"
#include "GlobalMotionRemover.h"
#include "CollectiveRandomPerturber.h"
#include "ForbidList.h"
#include <cmath>


Tools::Tools(
  RigidUnitSystem* sys_,
  ConstraintEnforcingPotential* cep_,
  NeighborTable* nt_,
  ProteinInfo* prot_ ) :
    sys(sys_),
    cep(cep_),
    nt(nt_),
    prot(prot_),
    revertReady( false ) {
  randomRotorPerturber = NULL;
  randomCenterPerturber = NULL;
  mom = NULL;
  momRandomRotorPerturber = NULL;
  momRandomCenterPerturber = NULL;
  globalMotionRemover = NULL;
  collectiveRandomPerturber = NULL;
}

Tools::~Tools() {
  delete randomRotorPerturber;
  delete randomCenterPerturber;
  delete mom;
  delete momRandomRotorPerturber;
  delete momRandomCenterPerturber;
  delete globalMotionRemover;
  delete collectiveRandomPerturber;
}

void Tools::setupRandomPert( double pertC, double pertR, double maxAngleRad ) {
  if ( randomRotorPerturber && randomCenterPerturber ) {
    randomCenterPerturber->setSize( pertC );
    randomRotorPerturber->setSize( pertR, maxAngleRad );
  }
  else {
    delete randomRotorPerturber;
    delete randomCenterPerturber;
    randomRotorPerturber = new RandomRotorPerturber( sys, pertR, maxAngleRad );
    randomCenterPerturber = new RandomCenterPerturber( sys, pertC );
  }
}

void Tools::setupMomentumPert( double pertC, double pertR, double momentumScaleFactor ) {
  delete mom;
  delete momRandomRotorPerturber;
  delete momRandomCenterPerturber;
  mom = new MomentumPerturber( sys );
  mom->setScaleFactor( momentumScaleFactor );
  momRandomRotorPerturber = new RandomRotorPerturber( sys, pertR );
  momRandomCenterPerturber = new RandomCenterPerturber( sys, pertC );
}

void Tools::setupCollectivePert() {
  collectiveRandomPerturber = new CollectiveRandomPerturber( prot, nt, sys, cep );
}

void Tools::collectivePert() {
  if ( collectiveRandomPerturber ) {
    cout << "  Doing Collective Perturbation" << endl;
    collectiveRandomPerturber->perturb();
  }
}


void Tools::randomPert() {
  if ( !randomCenterPerturber && !randomRotorPerturber ) return;

  cout << "  Doing random perturbation" << endl;
  if ( randomCenterPerturber ) randomCenterPerturber->perturb();
  if ( randomRotorPerturber ) randomRotorPerturber->perturb();
  minimize();
}

void Tools::randomPertRU( int ru ) {
  if ( randomCenterPerturber ) randomCenterPerturber->perturbRU( ru );
  if ( randomRotorPerturber ) randomRotorPerturber->perturbRU( ru );
}

void Tools::resetMomentumPert() {
  if ( mom ) mom->clear();
}

void Tools::momentumPert() {
  if ( mom && momRandomCenterPerturber && momRandomRotorPerturber ) {
    cout << "  Doing momentum perturbation" << endl;
    mom->determineDeltaQ();
    mom->setQ1();
    momRandomCenterPerturber->perturb();
    momRandomRotorPerturber->perturb();
    mom->perturb();
    minimize();
    if ( globalMotionRemover ) {
      globalMotionRemover->fitCurrentPointsToTarget();
    }
  }
}

void Tools::saveForRevert() {
  savedAbsolutePositions = sys->absolutePositions();
  revertReady = true;
}

void Tools::revertAll() {
  if ( !revertReady ) return;
  cout << "  Reverting structure." << endl;
  sys->setRigidUnits( savedAbsolutePositions );
  sys->update();
}

void Tools::revertProblemResidues( int& nPerturbedResi ) {
  //make a list of rigid units whose atoms are badly violating constraints
  nPerturbedResi = 0;
  if ( !cep || !prot ) return;

  set<int> problemRU;
  set<int> problemResi;
  double stretchthreshold = 0.15;

  if ( cep->bbhb ) {
    for ( BBHBContainer::const_iterator it = cep->bbhb->begin(); it != cep->bbhb->end(); it++ ) {
      if ( it->calcDistHA() - it->getConstraintMaxDistHA() > stretchthreshold ) {
        problemResi.insert( lookupResiFromAtom( it->geth() ) );
        problemResi.insert( lookupResiFromAtom( it->geta() ) );
      }
    }
  }

  if ( cep->hb ) {
    for ( HBContainer::const_iterator it = cep->hb->begin(); it != cep->hb->end(); it++ ) {
      if ( it->calcDistHA() - it->getConstraintMaxDistHA() > stretchthreshold ) {
        problemResi.insert( lookupResiFromAtom( it->geth() ) );
        problemResi.insert( lookupResiFromAtom( it->geta() ) );
      }
    }
  }

  if ( cep->ph ) {
    for ( PHContainer::const_iterator it = cep->ph->begin(); it != cep->ph->end(); it++ ) {
      if ( it->calcDist() - it->getCutoff() > stretchthreshold ) {
        problemResi.insert( lookupResiFromAtom( it->getp1() ) );
        problemResi.insert( lookupResiFromAtom( it->getp2() ) );
      }
    }
  }


  //rama
  double overlapthreshold = 0.07;
  if ( cep->rama ) {
    for ( RamaContainer::const_iterator it = cep->rama->begin(); it != cep->rama->end(); it++ ) {
      if ( it->violation() > overlapthreshold ) {
        problemResi.insert( lookupResiFromAtom( it->getp1() ) );
        problemResi.insert( lookupResiFromAtom( it->getp2() ) );
      }
    }
  }
  if ( cep->overridingRama ) {
    for ( MinDistConstraintContainer::const_iterator it = cep->overridingRama->begin();
          it != cep->overridingRama->end(); it++ ) {
      if ( it->violation() > overlapthreshold ) {
        problemResi.insert( lookupResiFromAtom( it->getp1() ) );
        problemResi.insert( lookupResiFromAtom( it->getp2() ) );
      }
    }
  }
  //side chain torsion
  if ( cep->sideChainTorsion ) {
    for ( SideChainTorsionContainer::const_iterator it = cep->sideChainTorsion->begin(); it != cep->sideChainTorsion->end(); it++ ) {
      if ( it->violation() > overlapthreshold ) {
        problemResi.insert( lookupResiFromAtom( it->getp1() ) );
        problemResi.insert( lookupResiFromAtom( it->getp2() ) );
      }
    }
  }
  if ( cep->overridingSC ) {
    for ( MinDistConstraintContainer::const_iterator it = cep->overridingSC->begin();
          it != cep->overridingSC->end(); it++ ) {
      if ( it->violation() > overlapthreshold ) {
        problemResi.insert( lookupResiFromAtom( it->getp1() ) );
        problemResi.insert( lookupResiFromAtom( it->getp2() ) );
      }
    }
  }
  //overlaps
  if ( cep->overlapEnergy ) {
    vector<MinDistConstraint>::const_iterator it = cep->overlapEnergy->constraints.begin();
    vector<MinDistConstraint>::const_iterator end = cep->overlapEnergy->constraints.end();
    for ( ; it != end; it++ ) {
      if ( it->violation() > overlapthreshold ) {
        problemResi.insert( lookupResiFromAtom( it->getp1() ) );
        problemResi.insert( lookupResiFromAtom( it->getp2() ) );
      }
    }
  }
  if ( cep->overridingMinDist ) {
    for ( MinDistConstraintContainer::const_iterator it = cep->overridingMinDist->begin();
          it != cep->overridingMinDist->end(); it++ ) {
      if ( it->violation() > overlapthreshold ) {
        problemResi.insert( lookupResiFromAtom( it->getp1() ) );
        problemResi.insert( lookupResiFromAtom( it->getp2() ) );
      }
    }
  }

  vector<int>::const_iterator begin;
  vector<int>::const_iterator end;
  vector<int>::const_iterator vecit;
  for ( set<int>::const_iterator resi = problemResi.begin(); resi != problemResi.end(); resi++ ) {
    for ( Resi::const_iterator atom = prot->resi( *resi ).begin();
          atom != prot->resi( *resi ).end(); atom++ ) {
      begin = sys->getRUlistFromP( atom->index() ).begin();
      end = sys->getRUlistFromP( atom->index() ).end();
      for ( vecit = begin; vecit != end; vecit++ ) {
        problemRU.insert( *vecit );
      }
    }
  }

  //revert each of the problem rigid units
  for ( set<int>::const_iterator ru = problemRU.begin(); ru != problemRU.end(); ru++ ) {
    sys->setRigidUnit( *ru, savedAbsolutePositions, sys->getRUPlistFromRU( *ru ) );
  }
  sys->collapseRotors();
  sys->update();

  /*
  //perturb each of the problem rigid units
  for ( set<int>::const_iterator ru = problemRU.begin(); ru != problemRU.end(); ru++ ) {
    randomPertRU( *ru );
  }
  sys->collapseRotors();
  sys->update();
  */
  nPerturbedResi = problemResi.size();
}

int Tools::lookupResiFromAtom( int a ) const {
  return prot->atom( a ).resi().index();
}

void Tools::setupGlobalMotionRemover_ReferenceInitialState() {
  globalMotionRemover = new GlobalMotionRemover( sys );

  //set the global fit mask
  vector<size_t> indexmask;
  size_t N = static_cast<size_t>( prot->natoms() );
  for ( size_t i = 0; i < N; i++ ) {
    //if ( prot->atom(i).name() == "CA" ) indexmask.push_back( i );
    indexmask.push_back( i );
  }
  globalMotionRemover->setIndexMask( indexmask );
  //set the current points as target, just once
  globalMotionRemover->setCurrentPointsAsTarget();

  //perform an initial fit
  globalMotionRemover->fitCurrentPointsToTarget();

}

void Tools::removeGlobalMotion() {
  if ( globalMotionRemover ) {
    globalMotionRemover->fitCurrentPointsToTarget();
  }
}

void Tools::problemPert() {
  int n = 0;

  revertProblemResidues( n );
  if ( n ) {
    cout << "  Reverted " << n << " problem residues" << endl;
    minimize();
  }
}

void Tools::minimize() {
  cep->minimize();
  cout << "  Performed " << cep->minim->getFinalStepNum() << " minimization cycles" << endl;
}

void Tools::outputSummaryLine() {
  cout << cep->generateSummaryString() << flush;
}

void TargetingTools::outputSummaryLine() {
  cout << ceptarg->generateSummaryString() << flush;
  int index;
  double dev;
  ceptarg->targetEnergy->calcMaxDeviation( index, dev );
  cout << "  Max Deviation " << index << " " << dev << endl;
}

void Tools::extraminimize() {
  cep->minim->setToleranceCondition( "maxPreconditionedGradComponent", 0.001 );
  minimize();
  cep->minim->setToleranceCondition( "maxPreconditionedGradComponent", 0.005 );
}

TargetingTools::TargetingTools(
  RigidUnitSystem* sys_,
  ConstraintEnforcingPotential_Targeting* cep_,
  NeighborTable* nt_,
  ProteinInfo* prot_ ) :
    Tools( sys_, cep_, nt_, prot_ ),
    ceptarg( cep_ ) {
  //constraintRemover_PQ = NULL;
  constraintRemover_Cutoff = NULL;
  constraintRemover_Random = NULL;
  rigidUnitReplace = NULL;
  forbidlist = NULL;
}

TargetingTools::~TargetingTools() {
  delete rigidUnitReplace;
  //delete constraintRemover_PQ;
  delete constraintRemover_Cutoff;
  delete constraintRemover_Random;
  delete forbidlist;
}

double TargetingTools::rmsd() {
  return ceptarg->targetEnergy->getRMSDtoTarget();
}
void TargetingTools::setRMSDconstraint( double c ) {
  ceptarg->targetEnergy->setRMSDconstraint( c );
}

void TargetingTools::adjustRMSDconstraint( double step ) {
  double initialrmsdConstraint = ceptarg->targetEnergy->getRMSDconstraint();
  double rmsdConstraint;
  rmsdConstraint = initialrmsdConstraint + step;
  if ( step < 0 && rmsdConstraint < 0 ) rmsdConstraint = 0;
  setRMSDconstraint( rmsdConstraint );
}

void TargetingTools::setupGlobalMotionRemover_Targeting() {
  globalMotionRemover = new GlobalMotionRemover( sys );
  vector<size_t> targetedIndices;
  vector<Vec3> targetPositions;
  ceptarg->targetEnergy->getTargetedIndices( targetedIndices );
  ceptarg->targetEnergy->getTargetPositions( targetPositions );
  globalMotionRemover->setIndexMask( targetedIndices );
  globalMotionRemover->setTarget( targetPositions );
  //perform an initial fit
  cout << "Performing fit to minimize RMSD of all targeted atoms" << endl;
  globalMotionRemover->fitCurrentPointsToTarget();
}

/*
void TargetingTools::setupConstraintRemover_PQ( double rmsd_upper, double rmsd_lower ) {
  constraintRemover_PQ = new ConstraintRemover_PQ( ceptarg );
  constraintRemover_PQ->initializePQ( rmsd_upper, rmsd_lower );
}

void TargetingTools::removeConstraints_PQ() {
  //remove any constraints that are scheduled for removal in the priority queue
  int nRemoved = 0;
  if ( constraintRemover_PQ ) {
    constraintRemover_PQ->remove( nRemoved );
    if ( nRemoved ) cout << "  Removed " << nRemoved << " constraints from queue" << endl;
  }
}
*/

void TargetingTools::setupReplace( double rmsd_upper, double rmsd_lower ) {
  if ( ceptarg->targneigh && ceptarg->targetEnergy && ceptarg->targmap ) {
    rigidUnitReplace = new RigidUnitReplace(
        sys, ceptarg->targetEnergy, &ceptarg->targmap->src2targ(), nt, ceptarg->targneigh );
    rigidUnitReplace->initializePQ( rmsd_upper, rmsd_lower );
  }
  else {
    cout << "Error, cannot initialize rigid unit replace," << endl;
    cout << "  because missing required target neighbor table or target energy or target map" << endl;
    exit(0);
  }
}

void TargetingTools::replace() {
  //rigid unit replace
  //Important: we must not do a replace when it has a chance
  //of being overwritten by a "revert" to saved coordinates

  if ( rigidUnitReplace ) {
    int count = 0;
    rigidUnitReplace->replace( count );
    if ( count ) {
      cout << "  Switched " << count << " rigid units to target geometry" << endl;
      minimize();
      outputSummaryLine();
      revertReady = false;
    }
  }

  /*
  //check for any problem shared points (do sparingly, possibly > 0.1 )
  //If any, replace all rigid units associated with these shared points and
  //reminimize.
  if ( rigidUnitReplace ) {
    double violationThreshold = 0.1;
    set<int> problemRU;
    set<int> additionalRU;
    set<int> additionalRU2;
    vector<int> violatingRUP;
    int replaced;
    do {
      problemRU.clear();
      additionalRU.clear();
      additionalRU2.clear();
      violatingRUP.clear();
      ceptarg->sharedPointsEnergy->getViolatingRigidUnitPoints( violationThreshold, violatingRUP );
      int n = violatingRUP.size();
      for ( int i = 0; i < n; i++ ) {
        int rup = violatingRUP[i];
        int p = sys->getPfromRUP( rup );
        problemRU.insert( sys->getRUlistFromP( p ).begin(), sys->getRUlistFromP( p ).end() );
      }
      if ( cep->overlapEnergy ) {
        vector<MinDistConstraint>::const_iterator it = cep->overlapEnergy->constraints.begin();
        vector<MinDistConstraint>::const_iterator end = cep->overlapEnergy->constraints.end();
        for ( ; it != end; it++ ) {
          if ( it->violation() > violationThreshold ) {
            problemRU.insert( sys->getRUlistFromP( it->getp1() ).begin(), sys->getRUlistFromP( it->getp1() ).end() );
            problemRU.insert( sys->getRUlistFromP( it->getp2() ).begin(), sys->getRUlistFromP( it->getp2() ).end() );
          }
        }
      }
      if ( cep->overridingMinDist ) {
        for ( MinDistConstraintContainer::const_iterator it = cep->overridingMinDist->begin();
              it != cep->overridingMinDist->end(); it++ ) {
          if ( it->violation() > violationThreshold ) {
            problemRU.insert( sys->getRUlistFromP( it->getp1() ).begin(), sys->getRUlistFromP( it->getp1() ).end() );
            problemRU.insert( sys->getRUlistFromP( it->getp2() ).begin(), sys->getRUlistFromP( it->getp2() ).end() );
          }
        }
      }

      for ( set<int>::const_iterator it_ru = problemRU.begin(); it_ru != problemRU.end(); it_ru++ ) {
        int ru = *it_ru;
        vector<int>::const_iterator it_p = sys->getPlistFromRU( ru ).begin();
        vector<int>::const_iterator p_end = sys->getPlistFromRU( ru ).end();
        for ( ; it_p != p_end; it_p++ ) {
          int p = *it_p;
          additionalRU.insert( sys->getRUlistFromP( p ).begin(), sys->getRUlistFromP( p ).end() );
        }
      }

      for ( set<int>::const_iterator it_ru = additionalRU.begin(); it_ru != additionalRU.end(); it_ru++ ) {
        int ru = *it_ru;
        vector<int>::const_iterator it_p = sys->getPlistFromRU( ru ).begin();
        vector<int>::const_iterator p_end = sys->getPlistFromRU( ru ).end();
        for ( ; it_p != p_end; it_p++ ) {
          int p = *it_p;
          additionalRU2.insert( sys->getRUlistFromP( p ).begin(), sys->getRUlistFromP( p ).end() );
        }
      }

      replaced = 0;
      for ( set<int>::const_iterator it_ru = additionalRU2.begin(); it_ru != additionalRU2.end(); it_ru++ ) {
        if ( rigidUnitReplace->isReplaceableRU( *it_ru ) ) {
          rigidUnitReplace->replaceRU( *it_ru );
          replaced++;
        }
      }

      if ( replaced ) {
        sys->update();
        cout << "  Switched " << replaced << " rigid units to target geometry" << endl;
        minimize();
        outputSummaryLine();
      }
    } while ( replaced > 0 );
  }
  */
}

void TargetingTools::addConstraints() {
  cout << "  Attempting to add constraints" << endl;
  if ( forbidlist ) tidyForbidList();
  if ( ceptarg->hbManager ) {
    ceptarg->hbManager->tighten();
    ceptarg->hbManager->findnew();
  }
  if ( ceptarg->phManager ) {
    ceptarg->phManager->tighten();
    ceptarg->phManager->findnew();
  }
}

void TargetingTools::setupConstraintRemover_Cutoff( 
    double stretchCutoff_bbhb, double stretchCutoff_hb, double stretchCutoff_ph ) {
  constraintRemover_Cutoff = new ConstraintRemover_Cutoff( ceptarg );
  constraintRemover_Cutoff->setStretchCutoff_bbhb( stretchCutoff_bbhb );
  constraintRemover_Cutoff->setStretchCutoff_hb( stretchCutoff_hb );
  constraintRemover_Cutoff->setStretchCutoff_ph( stretchCutoff_ph );
  cout << "Set up Cutoff Constraint Remover with stretch cutoff distances " << 
    stretchCutoff_bbhb << " (bbhb), " << stretchCutoff_hb << " (hb), " <<
    stretchCutoff_ph << " (ph)" << endl;
  if ( forbidlist ) constraintRemover_Cutoff->attachForbidList( forbidlist );

}

void TargetingTools::removeConstraints_Cutoff() {
  //check to see if any constraints are stretched
  if ( constraintRemover_Cutoff ) {
    int nRemoved;
    constraintRemover_Cutoff->remove( nRemoved );
    if ( nRemoved ) {
      cout << "  Removed " << nRemoved << " stretched constraints" << endl;
      //minimize();
    }
  }
}

void TargetingTools::setupConstraintRemover_Random( 
  double removeFrac_bbhb, double removeFrac_hb, double removeFrac_ph ) {
  constraintRemover_Random = new ConstraintRemover_Random( ceptarg );
  constraintRemover_Random->setRemoveFrac_bbhb( removeFrac_bbhb );
  constraintRemover_Random->setRemoveFrac_hb( removeFrac_hb );
  constraintRemover_Random->setRemoveFrac_ph( removeFrac_ph );
  cout << "Set up Random Constraint Remover with remove-frac " << 
    removeFrac_bbhb << " (bbhb), " << removeFrac_hb << " (hb), " <<
    removeFrac_ph << " (ph)" << endl;
  if ( forbidlist ) constraintRemover_Random->attachForbidList( forbidlist );
}

void TargetingTools::removeConstraints_Random() {
  //check to see if any constraints are stretched
  if ( constraintRemover_Random ) {
    int nRemoved;
    constraintRemover_Random->remove( nRemoved );
    if ( nRemoved ) {
      cout << "  Removed " << nRemoved << " randomly-chosen constraints" << endl;
      //minimize();
    }
  }
}

void TargetingTools::saveForRevert() {
  cout << "  Saving coordinates for Revert" << endl;
  savedAbsolutePositions = sys->absolutePositions();
  revertReady = true;
}

void TargetingTools::revertAll() {
  if ( !revertReady ) return;
  cout << "  Reverting structure." << endl;
  sys->setRigidUnits( savedAbsolutePositions );
  sys->update();
}

void TargetingTools::setupForbidList( int forbidTimeDuration_, Targeter* targeter_ ) {
  forbidlist = new ForbidList( forbidTimeDuration_, targeter_, prot, nt );
  if ( ceptarg->hbManager ) ceptarg->hbManager->attachForbidList( forbidlist );
  if ( ceptarg->phManager ) ceptarg->phManager->attachForbidList( forbidlist );
  if ( constraintRemover_Random ) constraintRemover_Random->attachForbidList( forbidlist );
  if ( constraintRemover_Cutoff ) constraintRemover_Cutoff->attachForbidList( forbidlist );
}

void TargetingTools::tidyForbidList() {
  if ( forbidlist ) forbidlist->tidy();
}

