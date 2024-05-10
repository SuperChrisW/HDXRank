#ifndef CONSTRAINTENFORCINGPOTENTIAL_H_
#define CONSTRAINTENFORCINGPOTENTIAL_H_

#include "SharedPoints.h"
#include "OverlapEnergy.h"
#include "SymmetryEnforcer.h"
#include "HBConstraint.h"
#include "HBManager.h"
#include "PHConstraint.h"
#include "PHManager.h"
#include "SideChainTorsion.h"
#include "Rama.h"
#include "SymmetryMatrices.h"
#include "Repulsion.h"
#include "MinimizeSystem.h"
#include "Energy.h"
#include "Gradient.h"
#include "Mismatch.h"
#include "GeneralizedCoords.h"
#include "MinDistConstraintContainer.h"
#include "Vec3.h"
class RigidUnitSystem;
class NeighborTable;
class Settings;
class ProteinInfo;
class PairTypeCutoffs;
#include <vector>

class ConstraintEnforcingPotential
{
public:
  ConstraintEnforcingPotential(
    RigidUnitSystem *rigidUnitSystem,
    const ProteinInfo& prot,
    const NeighborTable& nt,
    const Settings& settings );

  virtual ~ConstraintEnforcingPotential();

  void addEnergyTerm( EnergyTerm* term ) {
    energy_->addTerm( term );
  }
  void addGradientTerm( GradientTerm* term ) {
    gradient_->addTerm( term );
  }
  void addGradientTerm( GradientTerm_P* term ) {
    gradient_->addTerm( term );
  }
  void addMismatchTerm( MismatchTerm* term ) {
    mismatch_->addTerm( term );
  }

  double energy() { return energy_->calc(); }
  const GeneralizedCoords& gradient() { return gradient_->calc(); }
  const GeneralizedCoords& d2V_dQ2_diagonal() { return gradient_->calc_d2V_dQ2_diagonal(); }
  double maxMismatch() { return mismatch_->calc(); }

  void minimize() { minim->minimize(); }

  void reportConstraintViolations();
  void notifyTermsChanged() {
    energy_->notifyTermsChanged();
    gradient_->notifyTermsChanged();
  }

  void removeAllBreakableConstraints();

  void overrideProblemMinDistConstraints();
  void overrideConstraint( int p1, int p2, double newcutoff );

  virtual std::string generateColumnHeaderString();
  virtual std::string generateSummaryString();
  void writeConstraints();


  const std::vector<Vec3>& coords;
  SharedPoints *sharedPointsEnergy;
  OverlapEnergy *overlapEnergy;
  MinDistConstraintContainer *overridingMinDist;
  MinDistConstraintContainer *overridingSC;
  MinDistConstraintContainer *overridingRama;
  Repulsion *repulsion;
  PairTypeCutoffs *pairTypeCutoffs;
  BBHBContainer* bbhb;
  HBContainer* hb;
  HBManager* hbManager;
  PHContainer* ph;
  PHManager* phManager;
  SideChainTorsionContainer* sideChainTorsion;
  RamaContainer* rama;

  SymmetryEnforcer *symmetryEnforcerEnergy;
  SymmetryMatrices *symmetryMatrices;

  MinimizeSystem* minim;
private:
  Energy *energy_;
  Gradient *gradient_;
  Mismatch *mismatch_;

  void setupOverridingMinDist( std::string filename );

};

#endif /*CONSTRAINTENFORCINGPOTENTIAL_H_*/
