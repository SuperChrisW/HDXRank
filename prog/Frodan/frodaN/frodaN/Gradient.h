#ifndef GRADIENT_H_
#define GRADIENT_H_

#include <list>
#include <vector>
#include "GeneralizedCoords.h"
#include "Observable.h"
#include "Vec3.h"
class RigidUnitSystem;

class SecondDerivative {
public:
  double d2V_dx2;
  double d2V_dy2;
  double d2V_dz2;
  double d2V_dxdy;
  double d2V_dydz;
  double d2V_dzdx;
  SecondDerivative() {}
  SecondDerivative(
      double d2V_dx2, double d2V_dy2, double d2V_dz2,
      double d2V_dxdy, double d2V_dydz, double d2V_dzdx ) :
    d2V_dx2(d2V_dx2), d2V_dy2(d2V_dy2), d2V_dz2(d2V_dz2),
    d2V_dxdy(d2V_dxdy), d2V_dydz(d2V_dydz), d2V_dzdx(d2V_dzdx) {}
  SecondDerivative(const SecondDerivative& p) :
    d2V_dx2(p.d2V_dx2), d2V_dy2(p.d2V_dy2), d2V_dz2(p.d2V_dz2),
    d2V_dxdy(p.d2V_dxdy), d2V_dydz(p.d2V_dydz), d2V_dzdx(p.d2V_dzdx) {}
  ~SecondDerivative(){};

  SecondDerivative& operator=(const SecondDerivative& p) {
    d2V_dx2 = p.d2V_dx2; d2V_dy2 = p.d2V_dy2; d2V_dz2 = p.d2V_dz2;
    d2V_dxdy = p.d2V_dxdy; d2V_dydz = p.d2V_dydz; d2V_dzdx = p.d2V_dzdx;
    return *this;
  }

  void operator+=(const SecondDerivative& other) {
    d2V_dx2 += other.d2V_dx2; d2V_dy2 += other.d2V_dy2; d2V_dz2 += other.d2V_dz2;
    d2V_dxdy += other.d2V_dxdy; d2V_dydz += other.d2V_dydz; d2V_dzdx += other.d2V_dzdx;
  }
  void operator/=(double f) {
    d2V_dx2 /= f; d2V_dy2 /= f; d2V_dz2 /= f;
    d2V_dxdy /= f; d2V_dydz /= f; d2V_dzdx /= f;
  }
};

class Gradient;
class GradientTerm_P
{
public:
  GradientTerm_P() {}
  virtual ~GradientTerm_P() {}
  virtual void addToGradient_P( std::vector<Vec3> &dV_dr_Point,
      std::vector<SecondDerivative> &secondDerivative_Point ) = 0;
};

class GradientTerm
{
public:
  GradientTerm() {}
  virtual ~GradientTerm() {}
  virtual void addToGradient( std::vector<Vec3> &dV_dr_rigidUnitPoint,
      std::vector<SecondDerivative> &secondDerivative_rigidUnitPoint ) = 0;
};

class Gradient : public Observer {
public:
  Gradient( RigidUnitSystem *rigidUnitSystem_ );
  virtual ~Gradient();
  void addTerm( GradientTerm* term ) {
    gradientTerms.push_back(term);
    isUpdated = false;
    Observable *obs = dynamic_cast<Observable*>( term );
    if ( obs ) obs->registerObserver( this );
  }
  void addTerm( GradientTerm_P* term ) {
    gradientTerms_P.push_back(term);
    isUpdated = false;
    Observable *obs = dynamic_cast<Observable*>( term );
    if ( obs ) obs->registerObserver( this );
  }
  const GeneralizedCoords& calc() {
    if ( !isUpdated ) update();
    return gradientComponents;
  }
  const GeneralizedCoords& operator()() { return calc(); }
  const GeneralizedCoords& calc_d2V_dQ2_diagonal() {
    if ( !isUpdated ) update();
    return d2V_dQ2_diagonal;
  }
  void update();
  void receiveNotification( Observable *observable );
  void notifyTermsChanged() { isUpdated = false; }

private:
  RigidUnitSystem *rigidUnitSystem;
  std::list<GradientTerm*> gradientTerms;
  std::list<GradientTerm_P*> gradientTerms_P;

  //first derivative
  std::vector<Vec3> dV_dr_rigidUnitPoint;
  std::vector<Vec3> dV_dr_P;
  GeneralizedCoords gradientComponents;

  //second derivative
  std::vector<SecondDerivative> secondDerivative_rigidUnitPoint;
  std::vector<SecondDerivative> secondDerivative_P;
  GeneralizedCoords d2V_dQ2_diagonal;

  bool isUpdated;
  double rotormagcutoff2;
  void applyChainRule_P_to_RUP();
  void applyChainRule_RUP_to_RU();

};


#endif /*GRADIENT_H_*/
