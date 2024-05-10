#ifndef MOMENTUMPERTURBER_H_
#define MOMENTUMPERTURBER_H_

#include <vector>
class RigidUnitSystem;
#include "Vec3.h"
#include "Fit.h"

class MomentumPerturber
{
public:
  MomentumPerturber( RigidUnitSystem *rigidUnitSystem_ );
  virtual ~MomentumPerturber();
  void setScaleFactor( double s ) { scaleFactor = s; }
  void setQ1();
  void determineDeltaQ();
  void perturb();
  void clear();
  void activeMask( const std::vector<char>& activatedRU_ );

private:
  RigidUnitSystem *rigidUnitSystem;
  std::vector<Vec3> Q1_rigidUnitPoints;
  std::vector<char> activatedRU;
  std::vector<Vec3> step1_translation;
  std::vector<Vec3> step2_rotor;
  std::vector<Vec3> step2_centerOfRotation;
  Fit fit;
  bool isQ1set;
  bool readyToPerturb;
  double scaleFactor;
};

#endif /*MOMENTUMPERTURBER_H_*/
