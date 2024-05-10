#include "RandomRotorPerturber.h"
#include "RandomVector.h"
#include "RigidUnitSystem.h"
#include <algorithm>
#include "Vec3.h"
#include <cmath>
#include "mt19937ar.h"

using namespace std;

RandomRotorPerturber::RandomRotorPerturber( RigidUnitSystem *rigidUnitSystem_,
                                            double size_, double maxAngleRad_ ) :
  rigidUnitSystem( rigidUnitSystem_ ),
  size(size_),
  maxAngleRad( maxAngleRad_ ) {
  clear();
}

void RandomRotorPerturber::clear() {
  activatedRU.assign( rigidUnitSystem->nRigidUnits(), 1 );
}

void RandomRotorPerturber::activeMask( const vector<char>& activatedRU_ ) {
  if ( activatedRU_.size() != rigidUnitSystem->nRigidUnits() ) {
    cout << "Error RandomRotorPerturber: activeMask size does not match num. Rigid Units" << endl;
    exit(0);
  }
  activatedRU = activatedRU_;
}

void RandomRotorPerturber::perturbRU( int ru ) {
  if ( !activatedRU[ru] ) return;
  if ( rigidUnitSystem->hasZeroRadius(ru) ) return;

  Vec3 rotorPerturbation;
  double maxRotationAngle = min( maxAngleRad, size/rigidUnitSystem->radius(ru) );
  double randomRotationAngle = genrand_real1()*maxRotationAngle;
  double rotorMagnitude = 2.0 * sin(randomRotationAngle/2.0);
  generateRandomUnitVector( rotorPerturbation );
  rotorPerturbation *= rotorMagnitude;

  rigidUnitSystem->addToRotor( ru, rotorPerturbation );
}

void RandomRotorPerturber::perturb() {
  size_t nRU = rigidUnitSystem->nRigidUnits();

  for ( size_t ru = 0; ru < nRU; ru++ ) {
    perturbRU( ru );
  }
  rigidUnitSystem->collapseRotors();
  rigidUnitSystem->update();
}
