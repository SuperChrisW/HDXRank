#include "Gradient.h"
#include <cmath>
#include <iomanip>
#include <limits>
#include <cstdlib>
#include "Rotator.h"
#include "RigidUnitSystem.h"

using namespace std;

Gradient::Gradient( RigidUnitSystem *rigidUnitSystem_ ) :
  rigidUnitSystem( rigidUnitSystem_ ),
  isUpdated( false ),
  rotormagcutoff2( 1.0*1.0 )
{
}

Gradient::~Gradient()
{
}

void Gradient::receiveNotification( Observable *obs )
{
  isUpdated = false;
}


void Gradient::update() {
  //Initialize gradient to zero
  size_t nP = rigidUnitSystem->nPoints();
  size_t nRUP = rigidUnitSystem->nRigidUnitPoints();
  size_t nRU = rigidUnitSystem->nRigidUnits();

  dV_dr_P.assign( nP, Vec3( 0,0,0 ) );
  dV_dr_rigidUnitPoint.assign( nRUP, Vec3( 0,0,0 ) );
  gradientComponents.rotorsRU.assign( nRU, Vec3(0,0,0) );
  gradientComponents.centersRU.assign( nRU, Vec3(0,0,0) );

  //initialize second derivative to zero
  secondDerivative_P.assign( nP, SecondDerivative(0,0,0,0,0,0) );
  secondDerivative_rigidUnitPoint.assign( nRUP, SecondDerivative(0,0,0,0,0,0) );
  d2V_dQ2_diagonal.rotorsRU.assign( nRU, Vec3(0,0,0) );
  d2V_dQ2_diagonal.centersRU.assign( nRU, Vec3(0,0,0) );

  //if any rotor is big enough that it is in the unstable
  //region, collapse the rotor.
  for ( size_t ru=0; ru<nRU; ru++ ) {
    if ( rigidUnitSystem->rotors(ru).norm2() > rotormagcutoff2 ) {
      rigidUnitSystem->collapseRotor( ru );
    }
  }
  if ( rigidUnitSystem->BasePositionsChanged() ) {
    rigidUnitSystem->update();
  }

  //first do gradient calculation in terms of the cartesian positions of the atoms (mean points)
  for ( std::list<GradientTerm_P*>::iterator it = gradientTerms_P.begin();
        it != gradientTerms_P.end();
        it++ ) {
    // 'it' is an iterator to a pointer
    // *it is a pointer
    (*it)->addToGradient_P( dV_dr_P, secondDerivative_P );
  }

  //then apply chain rule to make with respect to cartesian positions of the rigid unit points
  if ( gradientTerms_P.size() != 0 ) {
    applyChainRule_P_to_RUP();
  }

  //add any gradient terms that are already in terms of the rigid unit points
  for ( std::list<GradientTerm*>::iterator it = gradientTerms.begin();
        it != gradientTerms.end();
        it++ ) {
     (*it)->addToGradient( dV_dr_rigidUnitPoint, secondDerivative_rigidUnitPoint );
  }

  //apply chain rule to make with respect to the rigid unit degrees of freedom
  applyChainRule_RUP_to_RU();

  isUpdated = true;
}

void Gradient::applyChainRule_P_to_RUP() {
  size_t nP = rigidUnitSystem->nPoints();
  for ( size_t p = 0; p < nP; p++ ) {
    const vector<int> *rup_list = &rigidUnitSystem->getRUPlistFromP(p);
    double N = static_cast<double>( rup_list->size() );
    Vec3 dV_dr = dV_dr_P[p];
    dV_dr /= N;
    SecondDerivative sd = secondDerivative_P[p];
    sd /= (N*N);
    int rup;
    for ( vector<int>::const_iterator rup_it = rup_list->begin();
          rup_it != rup_list->end();
          rup_it++ )
    {
      rup = *rup_it;

      dV_dr_rigidUnitPoint[rup] += dV_dr;
      secondDerivative_rigidUnitPoint[rup] += sd;
    }
  }
}

void Gradient::applyChainRule_RUP_to_RU() {
  size_t nRU = rigidUnitSystem->nRigidUnits();

  #pragma omp parallel
  {
    const Vec3 *dV_dr;
    const SecondDerivative *sD;
    Vec3 relativePosition;
    Rotator rotator;
    // since not all RUs have the same number of RUPs, it may be worthwhile to use
    // dynamic scheduling here
    #pragma omp for
    for ( size_t ru = 0; ru < nRU; ru++ ) {
      vector <int> rupList = rigidUnitSystem->getRUPlistFromRU(ru);
      for (size_t i = 0; i < rupList.size(); i++) {
	int rup = rupList[i];
        rotator.setRotor( rigidUnitSystem->rotors(ru) );
        dV_dr = &dV_dr_rigidUnitPoint[rup];
        //first derivative
        if ( ! ( abs(dV_dr_rigidUnitPoint[rup].x) < numeric_limits<double>::epsilon() &&
             abs(dV_dr_rigidUnitPoint[rup].y) < numeric_limits<double>::epsilon() &&
             abs(dV_dr_rigidUnitPoint[rup].z) < numeric_limits<double>::epsilon() ) )
        {

          gradientComponents.centersRU[ru] += *dV_dr;
          if ( !rigidUnitSystem->hasZeroRadius(ru) ) {

            relativePosition = rigidUnitSystem->absolutePositions(rup) -
                                rigidUnitSystem->centers(ru);
            // for rigid unit ru, the gradient of the rotor x-component is
            // the dV_dx*dx_dbx + dV_dy*dy_dbx + dV_dz*dz_dbx.
            gradientComponents.rotorsRU[ru].x +=
              dV_dr->x * rotator.dx_dbx( rigidUnitSystem->basePositions(rup) ) +
              dV_dr->y * rotator.dy_dbx( rigidUnitSystem->basePositions(rup) ) +
              dV_dr->z * rotator.dz_dbx( rigidUnitSystem->basePositions(rup) );

            // the gradient of the rotor y and z components follow the
            // same pattern
            gradientComponents.rotorsRU[ru].y +=
              dV_dr->x * rotator.dx_dby( rigidUnitSystem->basePositions(rup) ) +
              dV_dr->y * rotator.dy_dby( rigidUnitSystem->basePositions(rup) ) +
              dV_dr->z * rotator.dz_dby( rigidUnitSystem->basePositions(rup) );
            gradientComponents.rotorsRU[ru].z +=
              dV_dr->x * rotator.dx_dbz( rigidUnitSystem->basePositions(rup) ) +
              dV_dr->y * rotator.dy_dbz( rigidUnitSystem->basePositions(rup) ) +
              dV_dr->z * rotator.dz_dbz( rigidUnitSystem->basePositions(rup) );
          }
        }

        //second derivative
        if ( ! ( abs(secondDerivative_rigidUnitPoint[rup].d2V_dx2) < numeric_limits<double>::epsilon() &&
             abs(secondDerivative_rigidUnitPoint[rup].d2V_dy2) < numeric_limits<double>::epsilon() &&
             abs(secondDerivative_rigidUnitPoint[rup].d2V_dz2) < numeric_limits<double>::epsilon() &&
             abs(secondDerivative_rigidUnitPoint[rup].d2V_dxdy) < numeric_limits<double>::epsilon() &&
             abs(secondDerivative_rigidUnitPoint[rup].d2V_dydz) < numeric_limits<double>::epsilon() &&
             abs(secondDerivative_rigidUnitPoint[rup].d2V_dzdx) < numeric_limits<double>::epsilon() ) )
        {
          d2V_dQ2_diagonal.centersRU[ru].x += secondDerivative_rigidUnitPoint[rup].d2V_dx2;
          d2V_dQ2_diagonal.centersRU[ru].y += secondDerivative_rigidUnitPoint[rup].d2V_dy2;
          d2V_dQ2_diagonal.centersRU[ru].z += secondDerivative_rigidUnitPoint[rup].d2V_dz2;

          if ( !rigidUnitSystem->hasZeroRadius(ru) ) {

            relativePosition = rigidUnitSystem->absolutePositions(rup) -
                                rigidUnitSystem->centers(ru);
            sD = &secondDerivative_rigidUnitPoint[rup];
            d2V_dQ2_diagonal.rotorsRU[ru].x +=
              sD->d2V_dy2*relativePosition.z*relativePosition.z +
              sD->d2V_dz2*relativePosition.y*relativePosition.y -
              2.0*sD->d2V_dydz*relativePosition.y*relativePosition.z -
              dV_dr->y*relativePosition.y -
              dV_dr->z*relativePosition.z;
            d2V_dQ2_diagonal.rotorsRU[ru].y +=
              sD->d2V_dx2*relativePosition.z*relativePosition.z +
              sD->d2V_dz2*relativePosition.x*relativePosition.x -
              2.0*sD->d2V_dzdx*relativePosition.z*relativePosition.x -
              dV_dr->x*relativePosition.x -
              dV_dr->z*relativePosition.z;
            d2V_dQ2_diagonal.rotorsRU[ru].z +=
              sD->d2V_dx2*relativePosition.y*relativePosition.y +
              sD->d2V_dy2*relativePosition.x*relativePosition.x -
              2.0*sD->d2V_dxdy*relativePosition.y*relativePosition.x -
              dV_dr->y*relativePosition.y -
              dV_dr->x*relativePosition.x;
          }
        }
      }
    }
  } // end parallel region
}
