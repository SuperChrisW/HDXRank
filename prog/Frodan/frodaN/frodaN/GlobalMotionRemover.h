#ifndef GLOBALMOTIONREMOVER_H_
#define GLOBALMOTIONREMOVER_H_

using namespace std;

class RigidUnitSystem;
#include <vector>
#include "Vec3.h"
#include "Fit.h"

class GlobalMotionRemover
{
public:
  GlobalMotionRemover( RigidUnitSystem *rigidUnitSystem_ );
  virtual ~GlobalMotionRemover();
  void setIndexMask( const std::vector<size_t>& mask );
  void setTarget( const std::vector<Vec3>& target );
  void setCurrentPointsAsTarget();
  void fitCurrentPointsToTarget();
private:
  RigidUnitSystem *rigidUnitSystem;
  bool doMask;
  std::vector<size_t> mask;
  std::vector<Vec3> tempPoints;
  Fit fit;
};

#endif /*GLOBALMOTIONREMOVER_H_*/
