#include "RandomCenterPerturber.h"
#include "RandomVector.h"
#include "RigidUnitSystem.h"
#include "Vec3.h"
#include <cmath>
#include "mt19937ar.h"
#include <cstdlib>

using namespace std;

RandomCenterPerturber::RandomCenterPerturber( RigidUnitSystem *rigidUnitSystem_,
                                              double size_ ) :
  rigidUnitSystem( rigidUnitSystem_ ),
  size(size_) {
  clear();
}

void RandomCenterPerturber::clear() {
  activatedRU.assign( rigidUnitSystem->nRigidUnits(), 1 );
}

void RandomCenterPerturber::activeMask( const vector<char>& activatedRU_ ) {
  if ( activatedRU_.size() != rigidUnitSystem->nRigidUnits() ) {
    cout << "Error RandomCenterPerturber: activeMask size does not match num. Rigid Units" << endl;
    exit(0);
  }
  activatedRU = activatedRU_;
}

void RandomCenterPerturber::perturbRU( int ru ) {
  if ( !activatedRU[ru] ) return;
  Vec3 centerPerturbation;
  generateRandomUnitVector( centerPerturbation );
  centerPerturbation *= genrand_real2()*size;
  rigidUnitSystem->addToCenter( ru, centerPerturbation );
}

void RandomCenterPerturber::perturb() {
  size_t nRU = rigidUnitSystem->nRigidUnits();

  for ( size_t ru = 0; ru < nRU; ru++ ) {
    perturbRU( ru );
  }

  rigidUnitSystem->update();
}
