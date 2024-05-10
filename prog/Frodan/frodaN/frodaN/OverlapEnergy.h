/*
 * OverlapEnergy.h
 *
 *  Created on: Jul 24, 2009
 *      Author: dwfarrel
 */

#ifndef OVERLAPENERGY_H_
#define OVERLAPENERGY_H_

#include <vector>
#include "DistConstraint.h"
#include "Gradient.h"
#include "Vec3.h"
#include "Energy.h"
#include "Mismatch.h"
#include "VL.h"
#include "Repulsion.h"
class RigidUnitSystem;

class OverlapEnergy : public EnergyTerm,
                      public GradientTerm_P,
                      public MismatchTerm,
                      public Observer {
public:
  OverlapEnergy( RigidUnitSystem* sys, Repulsion* repulsion );
  virtual ~OverlapEnergy();

  double energy();

  void addToGradient_P(
      std::vector<Vec3> &dV_dr_Point,
      std::vector<SecondDerivative> &secondDerivative_Point );

  double mismatch();

  void clearPairs() { constraints.clear(); }
  void addPair( int p1, int p2, double cutoff ) {
    if ( constraints.size() > sizeTrigger ) {
      sizeTrigger = static_cast<size_t>( 1.1*constraints.size() );
      constraints.reserve( static_cast<size_t>( 1.2*constraints.size() ) );
    }
    constraints.push_back( MinDistConstraint( positions, k, p1, p2, cutoff ) );
  }

  void receiveNotification( Observable *obs ) {
    if ( sys == obs && !sys->AbsolutePositionsChanged() ) return;
    verlet.update();
  }

  void refresh() { verlet.makeNewList(); }

  double getk() const { return k; }

  std::vector<MinDistConstraint> constraints;


protected:
  RigidUnitSystem* sys;
  VL<Repulsion, OverlapEnergy> verlet;
  const std::vector<Vec3>* positions;
  double k;
  size_t sizeTrigger;
};

#endif /* OVERLAPENERGY_H_ */
