/*
 * Tools.h
 *
 *  Created on: Mar 24, 2010
 *      Author: dwfarrel
 */

#ifndef TOOLS_H_
#define TOOLS_H_

class RigidUnitSystem;
class RandomRotorPerturber;
class RandomCenterPerturber;
class MomentumPerturber;
class ConstraintEnforcingPotential;
class NeighborTable;
class ProteinInfo;
class GlobalMotionRemover;
class CollectiveRandomPerturber;

#include <vector>
#include "Vec3.h"

class Tools {
public:
  Tools(
    RigidUnitSystem* sys_,
    ConstraintEnforcingPotential* cep_,
    NeighborTable* nt_,
    ProteinInfo* prot_ );
  virtual ~Tools();

  void setupRandomPert( double pertC, double pertR, double maxAngleRad );
  void randomPert();
  void randomPertRU( int ru );

  void setupMomentumPert( double pertC, double pertR, double momentumScaleFactor = 1.0 );
  void resetMomentumPert();
  void momentumPert();

  void setupCollectivePert();
  void collectivePert();

  virtual void saveForRevert();
  virtual void revertAll();
  void problemPert();

  void setupGlobalMotionRemover_ReferenceInitialState();
  void removeGlobalMotion();

  void extraminimize();
  void minimize();

  virtual void outputSummaryLine();

protected:
  RigidUnitSystem* sys;
  ConstraintEnforcingPotential* cep;
  const NeighborTable* nt;
  const ProteinInfo* prot;
  RandomRotorPerturber* randomRotorPerturber;
  RandomCenterPerturber* randomCenterPerturber;
  MomentumPerturber* mom;
  RandomRotorPerturber* momRandomRotorPerturber;
  RandomCenterPerturber* momRandomCenterPerturber;
  GlobalMotionRemover* globalMotionRemover;
  CollectiveRandomPerturber* collectiveRandomPerturber;

  std::vector<Vec3> savedAbsolutePositions;
  bool revertReady;

  int lookupResiFromAtom( int a ) const;
  void revertProblemResidues( int& nPerturbedResi );

};

class ConstraintEnforcingPotential_Targeting;
class ConstraintRemover_Random;
class ConstraintRemover_Cutoff;
class RigidUnitReplace;
class ForbidList;
class Targeter;

class TargetingTools : public Tools {
public:
  TargetingTools(
    RigidUnitSystem* sys_,
    ConstraintEnforcingPotential_Targeting* cep_,
    NeighborTable* nt_,
    ProteinInfo* prot_ );
  virtual ~TargetingTools();

  virtual void saveForRevert();
  virtual void revertAll();

  void setupReplace( double rmsd_upper, double rmsd_lower );
  void replace();

  /*
  void setupConstraintRemover_PQ( double rmsd_upper, double rmsd_lower );
  void removeConstraints_PQ();
  */

  void setupConstraintRemover_Cutoff( 
    double stretchCutoff_bbhb, double stretchCutoff_hb, double stretchCutoff_ph );
  void removeConstraints_Cutoff();

  void setupConstraintRemover_Random( 
    double removeFrac_bbhb, double removeFrac_hb, double removeFrac_ph );
  void removeConstraints_Random();

  void setupForbidList( int forbidTimeDuration_, Targeter* targeter_ );
  void tidyForbidList();

  void setupGlobalMotionRemover_Targeting();

  void adjustRMSDconstraint( double step );
  double rmsd();
  void setRMSDconstraint( double c );

  //void attemptAddFinalConstraints();
  void addConstraints();

  virtual void outputSummaryLine();
protected:
  ConstraintEnforcingPotential_Targeting* ceptarg;
  //ConstraintRemover_PQ* constraintRemover_PQ;
  ConstraintRemover_Cutoff* constraintRemover_Cutoff;
  ConstraintRemover_Random* constraintRemover_Random;
  RigidUnitReplace* rigidUnitReplace;
  std::vector<char> isReplaceable;
  ForbidList* forbidlist;
};

#endif /* TOOLS_H_ */
