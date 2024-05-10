/*
 * Targeter.h
 *
 *  Created on: Aug 21, 2009
 *      Author: dwfarrel
 */

#ifndef TARGETER_H_
#define TARGETER_H_

class ProteinInfo;
class NeighborTable;
class RigidUnitSystem;
class ConstraintEnforcingPotential;
class ConstraintEnforcingPotential_Targeting;
class TargetingTools;
class Output;
class OutputFiles;
class Settings;
#include "DynamicConstraintRecorder.h"
#include "Vec3.h"

class Targeter {
public:

  Targeter( const Settings& settings );
  virtual ~Targeter();
  void run();
  double getRMSDtoTarget() const;
  int getIteration() const { return iter; }
private:
  ProteinInfo* prot;
  NeighborTable* nt;
  RigidUnitSystem* sys;
  ConstraintEnforcingPotential_Targeting* cep;
  TargetingTools* t;
  Output* output;
  OutputFiles* outputFiles;

  enum State { FORWARD, BACK, MOMENTUM, STOP };
  State state;

  size_t iter;
  size_t maxIter;
  int retryCount;
  int retryLimit;
  int randomPertFrequency;
  double delta0;
  double acceptableViolation;
  double violationTriggerExtraMinimization;
  bool doDynamicConstraints;
  int randomPertFreq;
  bool doReplace;
  bool doCollectivePert;
  bool doConstraintRemover_Cutoff;

  //backtracking variables
  bool doBacktrack;
  bool flag_noMoreBacktrack;
  double backtrackDelta;
  double backtrackDeltaInitial;
  double backtrackTarget;
  double backtrackDeltaLimit;
  double backtrackAcceptableMaxDev;
  double backtrackAcceptableRMSD;
  double backtrackDeltaScaleFactor;
  int momentumCount;
  int momentumLimit;
  double savedRMSD;

  DynamicConstraintRecorder<Targeter>* dynamicConstraintRecorder;

  bool structureAcceptable();
  void doExtraMinimizationIfNeeded();
  void outputStructure();

  //states
  void switchToForward();
  void forwardStep();
  void switchToBacktracking();
  void backwardsStep();
  void switchToMomentum();
  void momentumStep();
  void switchToStop();

  void initializeOptionalTargData(
    const Settings& settings,
    ProteinInfo *targprot,
    NeighborTable* targnt,
    RigidUnitSystem* targsys,
    ConstraintEnforcingPotential* targcep );
};

#endif /* TARGETER_H_ */
