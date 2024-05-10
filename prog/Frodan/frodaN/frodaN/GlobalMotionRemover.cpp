#include "GlobalMotionRemover.h"
#include "Fit.h"
#include "RigidUnitSystem.h"

void extractCoords(
    const vector<Vec3>& coords,
    const vector<size_t>& indices,
    vector<Vec3>& extractedCoords ) {
  const size_t n = indices.size();
  extractedCoords.resize( n );
  for ( size_t i = 0; i < n; i++ ) {
    extractedCoords[i] = coords[ indices[i] ];
  }
}

GlobalMotionRemover::GlobalMotionRemover(
      RigidUnitSystem *rigidUnitSystem_ ) :
  rigidUnitSystem( rigidUnitSystem_ ),
  doMask( false )
{
}

GlobalMotionRemover::~GlobalMotionRemover()
{
}

void GlobalMotionRemover::setIndexMask( const vector<size_t>& mask_ ) {
  mask = mask_;
  doMask = true;
}

void GlobalMotionRemover::setTarget( const vector<Vec3>& target ) {
  fit.setTargetAbsolutePoints( target );
}

void GlobalMotionRemover::setCurrentPointsAsTarget() {
  if ( doMask ) {
    extractCoords( rigidUnitSystem->meanPositions(), mask, tempPoints );
    fit.setTargetAbsolutePoints( tempPoints );
  }
  else {
    fit.setTargetAbsolutePoints( rigidUnitSystem->meanPositions() );
  }
}

void GlobalMotionRemover::fitCurrentPointsToTarget() {
  //set the source points
  if ( doMask ) {
    extractCoords( rigidUnitSystem->meanPositions(), mask, tempPoints );
    fit.setSourceAbsolutePoints( tempPoints );
  }
  else {
    fit.setSourceAbsolutePoints( rigidUnitSystem->meanPositions() );
  }

  //find the fit
  fit.simpleFit();
  Vec3 trans;
  Vec3 rotor;
  Vec3 centerOfRotation;
  fit.getFitStep1_translation( trans );
  fit.getFitStep2_rotor_centerOfRotation( rotor, centerOfRotation );

  //Apply the fit translation and rotation to the full set of coordinates.
  rigidUnitSystem->globalTranslate( trans );
  rigidUnitSystem->globalRotate( rotor, centerOfRotation );
  rigidUnitSystem->update();
}
