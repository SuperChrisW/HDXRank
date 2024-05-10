/*
 * DynamicConstraints.h
 *
 *  Created on: Jun 11, 2009
 *      Author: dwfarrel
 */

#ifndef DYNAMICCONSTRAINTS_H_
#define DYNAMICCONSTRAINTS_H_

class RigidUnitSystem;
class ConstraintEnforcingPotential;
class ProteinInfo;
class Settings;
class NeighborTable;
class Output;
class OutputFiles;
class Tools;
#include "Observable.h"
#include "HBManager.h"
#include "PHManager.h"
#include "HBConstraint.h"
#include "PHConstraint.h"

class DynamicConstraints {
public:
  DynamicConstraints( const Settings& settings );
  virtual ~DynamicConstraints();
  void run( int Ncycles );
  int getCycleCount() const { return cycle; }
private:
  RigidUnitSystem* sys;
  ConstraintEnforcingPotential* cep;
  ProteinInfo* prot;
  NeighborTable* nt;
  Output* output;
  OutputFiles* outputFiles;
  Tools *t;

  BBHBContainer* bbhb;
  HBContainer* hb;
  PHContainer* ph;
  HBManager* hbManager;
  PHManager* phManager;
  int cycle;

  bool dyn;
  bool doMomentumPert;
  bool doRandomPert;
  int switchToBreak;
  int switchToAdd;
  int switchOff;
  int switchBackToIterZero;

  void outputDynamicConstraintsColumnHeaders();
  void outputDynamicConstraintsColumns();
  void perturbSystemAndRestoreContraints();
  void addConstraints();

};

#endif /* DYNAMICCONSTRAINTS_H_ */
