/*
 * DistConstraint.h
 *
 *  Created on: Jul 24, 2009
 *      Author: dwfarrel
 */

#ifndef DISTCONSTRAINT_H_
#define DISTCONSTRAINT_H_

#include <vector>
#include "Vec3.h"
#include "Gradient.h"
#include "Energy.h"
#include <cmath>

class DistConstraint :
  public EnergyTerm,
  public GradientTerm_P {
public:
  DistConstraint();
  DistConstraint( const std::vector<Vec3>* positions_, double k_, int p1_, int p2_, double cutoff_ ) :
    positions(positions_),
    k(k_),
    p1(p1_),
    p2(p2_),
    cutoff( cutoff_ ) {}
  virtual ~DistConstraint();

  virtual double violation() const = 0;
  double calcDist2() const {
    return (*positions)[p1].dist2( (*positions)[p2] );
  }
  double calcDist() const { return sqrt( calcDist2() ); }
  double energy();
  void addToGradient_P(
    std::vector<Vec3> &dV_dr_P,
    std::vector<SecondDerivative> &secondDerivative_P );

  double getk() const { return k; }
  int getp1() const { return p1; }
  int getp2() const { return p2; }
  double getCutoff() const { return cutoff; }
  void setk( double k_ ) { k = k_; }
  void setPoints( const std::vector<Vec3>* positions_, int p1_, int p2_ ) {
    positions = positions_; p1 = p1_; p2 = p2_;
  }
  void setCutoff( double cutoff_ ) {
    cutoff = cutoff_;
  }

protected:
  virtual bool doesDistSquaredMeetConstraint( double dist2 ) const = 0;
  const std::vector<Vec3>* positions;
  double k;
  int p1;
  int p2;
  double cutoff;
};

class MinDistConstraint : public DistConstraint {
public:
  MinDistConstraint() {}
  MinDistConstraint( const std::vector<Vec3>* positions_, double k_, int p1_, int p2_, double cutoff_ ) :
    DistConstraint( positions_, k_, p1_, p2_, cutoff_ ) {}
  ~MinDistConstraint() {}
  double violation() const { return getCutoff() - calcDist(); }
private:
  virtual bool doesDistSquaredMeetConstraint( double dist2 ) const { return dist2 > cutoff*cutoff; }
};

class MaxDistConstraint : public DistConstraint {
public:
  MaxDistConstraint() {}
  MaxDistConstraint( const std::vector<Vec3>* positions_, double k_, int p1_, int p2_, double cutoff_ ) :
    DistConstraint( positions_, k_, p1_, p2_, cutoff_ ) {}
  ~MaxDistConstraint() {}
  double violation() const { return calcDist() - getCutoff(); }
private:
  virtual bool doesDistSquaredMeetConstraint( double dist2 ) const { return dist2 < cutoff*cutoff; }
};

#endif /* DISTCONSTRAINT_H_ */
