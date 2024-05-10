#include "MomentumPerturber.h"
#include "RigidUnitSystem.h"
#include <iostream>
#include <cstdlib>
#include <cmath>


using namespace std;

MomentumPerturber::MomentumPerturber( RigidUnitSystem *rigidUnitSystem_ ) :
  rigidUnitSystem( rigidUnitSystem_ ),
  scaleFactor( 1.0 )
{
  clear();
}

MomentumPerturber::~MomentumPerturber()
{
}

void MomentumPerturber::clear() {
  isQ1set = false;
  readyToPerturb = false;
  activatedRU.assign( rigidUnitSystem->nRigidUnits(), 1 );
}

void MomentumPerturber::activeMask( const vector<char>& activatedRU_ ) {
  if ( activatedRU_.size() != rigidUnitSystem->nRigidUnits() ) {
    cout << "Error MomentumPerturber: activeMask size does not match num. Rigid Units" << endl;
    exit(0);
  }
  activatedRU = activatedRU_;
}

void MomentumPerturber::setQ1() {
  Q1_rigidUnitPoints = rigidUnitSystem->absolutePositions();
  isQ1set = true;
}

void MomentumPerturber::determineDeltaQ() {
  if ( !isQ1set ) return;

  rigidUnitSystem->collapseRotors();
  size_t nRU = rigidUnitSystem->nRigidUnits();
  step1_translation.resize( nRU );
  step2_rotor.resize( nRU );
  step2_centerOfRotation.resize( nRU );
  for ( size_t ru = 0; ru < nRU; ru++ ) {
    vector<int> ruplist;
    if ( !activatedRU[ru] ) continue;
    ruplist = rigidUnitSystem->getRUPlistFromRU( ru );
    fit.setSourceAbsolutePoints( Q1_rigidUnitPoints, ruplist );
    fit.setTargetAbsolutePoints( rigidUnitSystem->absolutePositions(), ruplist );
    fit.simpleFit();
    fit.getFitStep1_translation( step1_translation[ru] );
    fit.getFitStep2_rotor_centerOfRotation( step2_rotor[ru], step2_centerOfRotation[ru] );
  }

  isQ1set = false;
  readyToPerturb = true;
}

void MomentumPerturber::perturb() {
  if ( !readyToPerturb ) return;

  size_t nRU = rigidUnitSystem->nRigidUnits();
  const double maglim = 0.5;
  const double maglim2 = maglim*maglim;
  double mag2;

  rigidUnitSystem->collapseRotors();

  for ( size_t ru = 0; ru < nRU; ru++ ) {
    if ( !activatedRU[ru] ) continue;

    //get the net translation of this rigid unit from state Q1 to state Q2
    Vec3 translation = step1_translation[ru]*scaleFactor;

    //if translation was too large, truncate
    mag2 = translation.norm2();
    if ( mag2 > maglim2 ) {
      translation *= maglim/sqrt(mag2);
    }

    //apply perturbation
    rigidUnitSystem->addToCenter( ru, translation );

    if ( rigidUnitSystem->hasZeroRadius(ru) ) continue;

    //get the net rotation of this rigid unit from state Q1 to Q2
    Vec3 rotor = step2_rotor[ru]*scaleFactor; //note this scaling is not perfect because we really should be scaling the rotation angle rather than the rotor magnitude.
    Vec3 centerOfRotation = step2_centerOfRotation[ru];

    //if rotation was too large, truncate.
    //use the rotor magnitude (which is the rotation angle for small angles),
    //scaled by the radius
    mag2 = rotor.norm2()*rigidUnitSystem->radius(ru)*rigidUnitSystem->radius(ru);
    if ( mag2 > maglim2 ) {
      rotor *= maglim/sqrt(mag2);
    }

    // apply perturbation
    Rotator rotator;
    rotator.setRotor( rotor );
    rigidUnitSystem->rotate( ru, rotator, centerOfRotation );
  }
  rigidUnitSystem->collapseRotors();
  rigidUnitSystem->update();

  readyToPerturb = false;
}
