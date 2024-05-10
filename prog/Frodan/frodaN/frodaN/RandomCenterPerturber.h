#ifndef RANDOMCENTERPERTURBER_H_
#define RANDOMCENTERPERTURBER_H_

class RigidUnitSystem;
#include <vector>

class RandomCenterPerturber
{
public:
  RandomCenterPerturber( RigidUnitSystem *rigidUnitSystem_,
                         double size_ );
	virtual ~RandomCenterPerturber() {}

  void perturb();
  void perturbRU( int ru );
  void clear();
  void activeMask( const std::vector<char>& activatedRU_ );
  void setSize( double size_ ) { size = size_; }

private:
  RigidUnitSystem *rigidUnitSystem;
  std::vector<char> activatedRU;
  double size;
};

#endif /*RANDOMCENTERPERTURBER_H_*/
