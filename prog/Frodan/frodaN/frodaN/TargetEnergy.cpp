#include "TargetEnergy.h"
#include "RigidUnitSystem.h"
#include <cmath>
#include <limits>

using namespace std;

TargetEnergy::TargetEnergy( RigidUnitSystem *rigidUnitSystem_ ) :
  enabled( true ),
  rmsdConstraint(0),
  C(0),
  sqrtSumSquares(0),
  rigidUnitSystem(rigidUnitSystem_),
  k( 10 ),
  forwards( true )
{
  rigidUnitSystem->registerObserver( this );
  mapPtoTargetIndex.resize( rigidUnitSystem->nPoints(), -1 );
}

TargetEnergy::~TargetEnergy()
{
}

void TargetEnergy::update() {
  //calculate difference vector from target point to source (mean) point
  sqrtSumSquares = 0;
  for ( size_t i = 0; i < targets.size(); i++ ) {
    targets[i].differenceVec = rigidUnitSystem->meanPositions( targets[i].p );
    targets[i].differenceVec -= targets[i].targetVec;
    sqrtSumSquares += targets[i].differenceVec.norm2();
  }
  sqrtSumSquares = sqrt( sqrtSumSquares );
  notifyObservers();
}

void TargetEnergy::calcMaxDeviation( int& atomindex, double& dev ) const {
  size_t N = targets.size();
  double maxDev2 = 0;
  int maxi = 0;
  for ( size_t i = 0; i < N; i++ ) {
    double dev2 = targets[i].differenceVec.norm2();
    if ( dev2 > maxDev2 ) {
      maxDev2 = dev2;
      maxi = i;
    }
  }

  atomindex = targets[maxi].p;
  dev = sqrt( maxDev2 );
}

void TargetEnergy::setDirectionOfReactionCoordinate() {
  size_t N = targets.size();
  u.resize( N );
  double norm2 = 0;
  for ( size_t i = 0; i < N; i++ ) {
    u[i].x = -targets[i].differenceVec.x;
    u[i].y = -targets[i].differenceVec.y;
    u[i].z = -targets[i].differenceVec.z;
    norm2 += u[i].norm2();
  }
  double norm = sqrt(norm2);
  for ( size_t i = 0; i < N; i++ ) {
    u[i] /= norm;
  }
}

void TargetEnergy::getProjections( double &coord1, double &coord2, double &rmsd ) const {
  size_t N = targets.size();
  coord1 = 0;
  for ( size_t i = 0; i < N; i++ ) {
    coord1 += u[i].dot( targets[i].differenceVec );
  }
  coord1 /= sqrt( static_cast<double>(N) );
  rmsd = getRMSDtoTarget();
  coord2 = sqrt( rmsd*rmsd - coord1*coord1 );
}

void TargetEnergy::getTargetedIndices( std::vector<size_t>& indices ) const {
  size_t N = targets.size();
  indices.resize( N );
  for ( size_t i = 0; i < N; i++ ) {
    indices[i] = targets[i].p;
  }
}
void TargetEnergy::getTargetPositions( std::vector<Vec3>& targetPositions ) const {
  size_t N = targets.size();
  targetPositions.resize( N );
  for ( size_t i = 0; i < N; i++ ) {
    targetPositions[i] = targets[i].targetVec;
  }
}


void TargetEnergy::addToGradient( std::vector<Vec3> &dV_dr_rigidUnitPoint,
                                  std::vector<SecondDerivative> &secondDerivative_rigidUnitPoint ) {
  if ( !enabled ||
       forwards && sqrtSumSquares <= C ||
       !forwards && sqrtSumSquares >= C ) return;

  double sumSquares = sqrtSumSquares*sqrtSumSquares;
  SecondDerivative secondDerivative;

  bool CequalsZero = rmsdConstraint <= numeric_limits<double>::epsilon();

  //The 1st and 2nd derivatives take two analytical forms depending on whether C is zero.
  //The C equals zero case is equivalent to harmonic restraints from
  //the atoms to their target positions, and the derivatives take on a simpler form.
  //The general case (C not necessarily zero) corresponds to a harmonic
  //restraint on the value of the overall RMSD to the target, and the derivatives
  //are a little more complicated.

  //First, the C equals zero case.
  if ( CequalsZero ) {
    for ( size_t i = 0; i < targets.size(); i++ ) {
      const vector<int> *rup_list;
      rup_list = &rigidUnitSystem->getRUPlistFromP( targets[i].p );
      Vec3 grad = targets[i].differenceVec;
      grad *= k/static_cast<double>( rup_list->size() );
      secondDerivative.d2V_dx2 = secondDerivative.d2V_dy2 =
        secondDerivative.d2V_dz2 = k/static_cast<double>( rup_list->size()*rup_list->size() );

      int rup;
      for ( vector<int>::const_iterator rup_it = rup_list->begin();
            rup_it != rup_list->end();
            rup_it++ )
      {
        rup = *rup_it;
        dV_dr_rigidUnitPoint[rup] += grad;
        secondDerivative_rigidUnitPoint[rup].d2V_dx2 += secondDerivative.d2V_dx2;
        secondDerivative_rigidUnitPoint[rup].d2V_dy2 += secondDerivative.d2V_dy2;
        secondDerivative_rigidUnitPoint[rup].d2V_dz2 += secondDerivative.d2V_dz2;
        //the cross terms, d2V_dxdy etc., are zero in this case, so no need to
        //add them to the second derivative
      }
    }
    return;
  }

  //Otherwise, do the general case (C nonzero)
  double C_over_sqrtSumSquares = C/sqrtSumSquares;
  double factor1a = k*(1.0 - C_over_sqrtSumSquares);
  for ( size_t i = 0; i < targets.size(); i++ ) {
    const vector<int> *rup_list;
    rup_list = &rigidUnitSystem->getRUPlistFromP( targets[i].p );
    Vec3 grad = targets[i].differenceVec;
    grad *= factor1a/static_cast<double>( rup_list->size() );

    double factor2a = k/static_cast<double>( rup_list->size()*rup_list->size() );
    secondDerivative.d2V_dx2 = factor2a*( 1.0 - C_over_sqrtSumSquares*
        ( 1.0 - targets[i].differenceVec.x*targets[i].differenceVec.x/sumSquares ) );
    secondDerivative.d2V_dy2 = factor2a*( 1.0 - C_over_sqrtSumSquares*
        ( 1.0 - targets[i].differenceVec.y*targets[i].differenceVec.y/sumSquares ) );
    secondDerivative.d2V_dz2 = factor2a*( 1.0 - C_over_sqrtSumSquares*
        ( 1.0 - targets[i].differenceVec.z*targets[i].differenceVec.z/sumSquares ) );
    secondDerivative.d2V_dxdy = factor2a*C_over_sqrtSumSquares*
        targets[i].differenceVec.x*targets[i].differenceVec.y/sumSquares;
    secondDerivative.d2V_dydz = factor2a*C_over_sqrtSumSquares*
        targets[i].differenceVec.y*targets[i].differenceVec.z/sumSquares;
    secondDerivative.d2V_dzdx = factor2a*C_over_sqrtSumSquares*
        targets[i].differenceVec.z*targets[i].differenceVec.x/sumSquares;

    int rup;
    for ( vector<int>::const_iterator rup_it = rup_list->begin();
          rup_it != rup_list->end();
          rup_it++ )
    {
      rup = *rup_it;
      dV_dr_rigidUnitPoint[rup] += grad;
      secondDerivative_rigidUnitPoint[rup].d2V_dx2 += secondDerivative.d2V_dx2;
      secondDerivative_rigidUnitPoint[rup].d2V_dy2 += secondDerivative.d2V_dy2;
      secondDerivative_rigidUnitPoint[rup].d2V_dz2 += secondDerivative.d2V_dz2;
      secondDerivative_rigidUnitPoint[rup].d2V_dxdy += secondDerivative.d2V_dxdy;
      secondDerivative_rigidUnitPoint[rup].d2V_dydz += secondDerivative.d2V_dydz;
      secondDerivative_rigidUnitPoint[rup].d2V_dzdx += secondDerivative.d2V_dzdx;
    }
  }
}

