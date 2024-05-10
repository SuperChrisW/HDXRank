#ifndef RANDOMROTORPERTURBER_H_
#define RANDOMROTORPERTURBER_H_

class RigidUnitSystem;
#include <vector>
#include <cmath>

class RandomRotorPerturber
{
public:
  RandomRotorPerturber( RigidUnitSystem *rigidUnitSystem_,
                        double size_, double maxAngleRad_ = M_PI );
  virtual ~RandomRotorPerturber() {}

  void perturb();
  void perturbRU( int ru );
  void clear();
  void activeMask( const std::vector<char>& activatedRU_ );
  void setSize( double size_, double maxAngleRad_ = M_PI ) {
    size = size_;
    maxAngleRad = maxAngleRad_;
  }

private:
  RigidUnitSystem *rigidUnitSystem;
  std::vector<char> activatedRU;
  double size;
  double maxAngleRad;
};

#endif /*RANDOMROTORPERTURBER_H_*/
