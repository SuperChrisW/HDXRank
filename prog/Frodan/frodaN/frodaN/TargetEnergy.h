#ifndef TARGETENERGY_H_
#define TARGETENERGY_H_

class RigidUnitSystem;
#include "Vec3.h"
#include "Energy.h"
#include "Gradient.h"
#include "Observable.h"
#include "RigidUnitSystem.h"
#include <vector>
#include <cmath>

class TargetData {
public:
  TargetData() {}
  TargetData( int p_, const Vec3& targetVec_ ) {
    p = p_;
    targetVec = targetVec_;
  }
  ~TargetData() {}

  int p;
  Vec3 targetVec;
  Vec3 differenceVec;
};

class TargetEnergy : public EnergyTerm,
                     public GradientTerm,
                     public Observer,
                     public Observable
{
public:
	TargetEnergy( RigidUnitSystem *rigidUnitSystem );

	virtual ~TargetEnergy();
	void addTargetPoint( int p, const Vec3& targetVec ) {
    mapPtoTargetIndex[p] = targets.size();
	  targets.push_back( TargetData( p, targetVec ) );
	}
	void clear() {
	  targets.clear();
	  mapPtoTargetIndex.clear();
	  u.clear();
    update();
	}
  void setDirectionOfReactionCoordinate();

  void enable() {
    enabled = true;
    update();
  }

  bool isEnabled() const { return enabled; }

  void disable() {
    enabled = false;
    notifyObservers();
  }

  double energy() {
    if ( !enabled ) return 0;
    if ( forwards && sqrtSumSquares > C || !forwards && sqrtSumSquares < C )
      return 0.5*k*(sqrtSumSquares - C)*(sqrtSumSquares - C);
    else
      return 0;
  }
  void addToGradient( std::vector<Vec3> &dV_dr_rigidUnitPoint,
                      std::vector<SecondDerivative> &secondDerivative_rigidUnitPoint );
  void receiveNotification( Observable *obs ) {
    RigidUnitSystem *sys = dynamic_cast<RigidUnitSystem*>( obs );
    if ( sys && sys->AbsolutePositionsChanged() ) update();
    else if ( !sys ) update();
  }
  void setRMSDconstraint( double rmsdConstraint_ ) {
    rmsdConstraint = rmsdConstraint_;
    C = rmsdConstraint*sqrt( static_cast<double>(targets.size()) );
    enable();
  }
  double getRMSDconstraint() const {
    return rmsdConstraint;
  }
  double getRMSDtoTarget() const {
    return sqrtSumSquares/sqrt( static_cast<double>(targets.size()) );
  }

  double calcMaxDeviation() const {
    int index;
    double dev;
    calcMaxDeviation( index, dev );
    return dev;
  }
  void calcMaxDeviation( int& atomindex, double& dev ) const;
  void getProjections( double &coord1, double &coord2, double &rmsd ) const;
  void getTargetedIndices( std::vector<size_t>& indices ) const;
  void getTargetPositions( std::vector<Vec3>& targetPositions ) const;
  void getTargetPosition( int p, bool& found, Vec3& position ) const {
    if ( mapPtoTargetIndex[p] == -1 ) {
      found = false;
    }
    else {
      found = true;
      position = targets[ mapPtoTargetIndex[p] ].targetVec;
    }
  }

  double getDistToTarget( int p ) const {
    Vec3 position;
    bool found;
    getTargetPosition( p, found, position );
    return found ? sqrt( rigidUnitSystem->meanPositions(p).dist2( position ) ) : 0;
  }

  //updates the vector differences from source points to target points,
  //and updates the value of sqrtSumSquares given the current targets
  void update();

  void setForwards() {
    forwards = true;
    notifyObservers();
  }
  void setBackwards() {
    forwards = false;
    notifyObservers();
  }

  void getDistSquared( int p, bool& found, double& dist2 ) const {
    if ( mapPtoTargetIndex[p] == -1 ) {
      found = false;
      dist2 = 0;
    }
    else {
      found = true;
      dist2 = targets[ mapPtoTargetIndex[p] ].differenceVec.norm2();
    }
  }
private:
  bool enabled;
  double rmsdConstraint;
  double C;
  double sqrtSumSquares;
  RigidUnitSystem *rigidUnitSystem;
  std::vector<TargetData> targets;
  double k;
  std::vector<Vec3> u;
  std::vector<int> mapPtoTargetIndex;
  bool forwards;
};

#endif /*TARGETENERGY_H_*/
