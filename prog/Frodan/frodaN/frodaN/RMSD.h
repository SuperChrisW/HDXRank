#ifndef RMSD_H_
#define RMSD_H_

#include "Vec3.h"
#include "RigidUnitSystem.h"
#include <vector>

class RMSD
{
public:
	RMSD( const RigidUnitSystem *rigidUnitSystem_ ) :
	  rigidUnitSystem( rigidUnitSystem_ ),
	  referencePoints( rigidUnitSystem->meanPositions() ) {}

	virtual ~RMSD() {}

	void setCurrentPointsAsReference() {
	  referencePoints = rigidUnitSystem->meanPositions();
	}

	double calcRMSD() {
	  return calcRMSD( referencePoints, rigidUnitSystem->meanPositions() );
	}

	double calcRMSD( const std::vector<Vec3>& points1, const std::vector<Vec3>& points2 ) const;

private:
  const RigidUnitSystem *rigidUnitSystem;
  std::vector<Vec3> referencePoints;
};

#endif /*RMSD_H_*/
