/*
 * Swap.cpp
 *
 *  Created on: May 15, 2009
 *      Author: dwfarrel
 */

#include "Swap.h"
#include <cmath>
#include "Fit.h"
#include "Rotator.h"
#include <iostream>

using namespace std;

SwapBase::~SwapBase() {}

SwapSet_2Fold::SwapSet_2Fold() : SwapBase() {
  altmaps.resize( 1 );
}

void SwapSet_2Fold::clear() {
  indices.clear();
  altmaps[0].clear();
}

void SwapSet_2Fold::addPair( int i, int j ) {
  indices.push_back( i );
  indices.push_back( j );

  //altmaps[0] gets the following map:
  // [0] -> 1
  // [1] -> 0
  //-----------
  // [2] -> 3
  // [3] -> 2
  //-----------
  // ...etc.
  int n = altmaps[0].size();
  altmaps[0].push_back( n+1 );
  altmaps[0].push_back( n );
}

SwapPoint_3Fold::SwapPoint_3Fold() : SwapBase() {
  altmaps.resize( 2 );

  altmaps[0].resize( 3 );
  altmaps[0][0] = 1;
  altmaps[0][1] = 2;
  altmaps[0][2] = 0;

  altmaps[1].resize( 3 );
  altmaps[1][0] = 2;
  altmaps[1][1] = 0;
  altmaps[1][2] = 1;
}

void SwapPoint_3Fold::setTriple( int i, int j, int k ) {
  indices.resize( 3 );
  indices[0] = i;
  indices[1] = j;
  indices[2] = k;
}

void SwapPoint_3Fold::setTriple( const std::vector<int>& triple ) {
  indices = triple;
  indices.resize(3);
}

CoordinateSwapper::CoordinateSwapper(
    vector<Vec3>& state1coords_,
    const vector<Vec3>& state2coords_,
    const map<int,int>* atommap_ ) :
    state1coords( state1coords_ ),
    state2coords( state2coords_ ),
    atommap( atommap_ )
  {
  }

void CoordinateSwapper::swap(
    const vector<int>& fittingSubset1,
    const SwapBase& swapBase ) {

    const vector<int>& swappingSubset1 = swapBase.getIndices();

    //Get the converted index lists for structure 2.
    //If any of the structure 1 atoms do not have counterparts in structure 2,
    //then we do not proceed.
    bool success = false;
    vector<int> fittingSubset2;
    conv1to2( fittingSubset1, fittingSubset2, success );
    if ( !success ) return;
    vector<int> swappingSubset2;
    conv1to2( swappingSubset1, swappingSubset2, success );
    if ( !success ) return;

    //extract coordinates in structures 1 and 2 for the "fitting subset",
    //and find the best fit that takes the structure 2 coords to the structure 1 coords
    vector<Vec3> coordsSubset1;
    vector<Vec3> coordsSubset2;
    extractCoords( state1coords, fittingSubset1, coordsSubset1 );
    extractCoords( state2coords, fittingSubset2, coordsSubset2 );
    fit.setSourceAbsolutePoints( coordsSubset2 );
    fit.setTargetAbsolutePoints( coordsSubset1 );
    fit.simpleFit();
    Vec3 trans;
    Vec3 rotor;
    Vec3 centerOfRotation;
    Rotator rotator;
    fit.getFitStep1_translation( trans );
    fit.getFitStep2_rotor_centerOfRotation( rotor, centerOfRotation );
    rotator.setRotor( rotor );

    //Now extract the coordinates in structures 1 and 2 for the "swapping subset".
    //Apply the fit translation and rotation to the "swapping subset".
    extractCoords( state1coords, swappingSubset1, coordsSubset1 );
    extractCoords( state2coords, swappingSubset2, coordsSubset2 );
    const int n = coordsSubset2.size();
    for ( int i = 0; i < n; i++ ) {
      coordsSubset2[i] += trans;
      rotator.rotateAboutCenter( coordsSubset2[i], centerOfRotation, coordsSubset2[i] );
    }

    //now find the best swap
    const vector< vector<int> >& altmaps = swapBase.getAltMaps();
    pickBest( swappingSubset1, altmaps, coordsSubset2 );
}

void CoordinateSwapper::conv1to2(
    const vector<int>& indices1,
    vector<int>& indices2,
    bool& success ) const {

    success = true;

    if ( !atommap ) {
      indices2 = indices1;
      return;
    }

    indices2.resize( indices1.size() );

    size_t n = indices1.size();
    map< int, int >::const_iterator it;
    for ( size_t i = 0; i < n; i++ ) {
      it = atommap->find( indices1[i] );
      if ( it != atommap->end() ) {
        indices2[i] = it->second;
      }
      else {
        success = false;
        indices2.clear();
        return;
      }
    }
}

void CoordinateSwapper::swapIndices(
    const vector<int>& indices,
    const vector<int>& remapping,
    vector<int>& swappedIndices ) const {
    int n = indices.size();
    if ( indices.size() != remapping.size() ) {
      cout << "Error, sizes do not match" << endl;
      exit(0);
    }
    swappedIndices.resize( indices.size() );
    for ( int i = 0; i < n; i++ ) {
      swappedIndices[i] = indices[ remapping[i] ];
    }
}

void CoordinateSwapper::extractCoords(
    const vector<Vec3>& coords,
    const vector<int>& indices,
    vector<Vec3>& extractedCoords ) const {
    int n = indices.size();
    extractedCoords.resize( indices.size() );
    for ( int i = 0; i < n; i++ ) {
      extractedCoords[i] = coords[ indices[i] ];
    }
}

void CoordinateSwapper::setCoords(
    vector<Vec3>& coords,
    const vector<int>& indices,
    const vector<Vec3>& newCoords ) const {
    int n = indices.size();
    for ( int i = 0; i < n; i++ ) {
      coords[ indices[i] ] = newCoords[i];
    }
}

double CoordinateSwapper::energy( const vector<Vec3>& positions1, const vector<Vec3>& positions2 ) const {
    double E = 0;
    size_t N = positions1.size();
    if ( positions2.size() != N ) {
      cout << "Error in coordinate swapper: vectors different sizes" << endl;
      exit(0);
    }
    for ( size_t i = 0; i < N; i++ ) {
      E += positions1[i].dist2( positions2[i] );
    }
    return E;
}

void CoordinateSwapper::pickBest( const vector<int>& indices,
    const vector< vector<int> >& altmaps,
    const vector<Vec3>& reference ) {

    //try out all the alternative mappings, see which matches
    //best with the reference positions
    double minE = numeric_limits<double>::max();
    int bestSwapID = -1;
    vector<int> swappedIndices;
    vector<Vec3> swappedCoords;
    const int nAlternates = altmaps.size();
    for ( int i = 0; i < nAlternates; i++ ) {
      swapIndices( indices, altmaps[i], swappedIndices );
      extractCoords( state1coords, swappedIndices, swappedCoords );
      double E = energy( swappedCoords, reference );
      if ( E < minE ) {
        minE = E;
        bestSwapID = i;
      }
    }
    //still must try out the original mapping
    extractCoords( state1coords, indices, swappedCoords );
    double Eorig = energy( swappedCoords, reference );

    //If the original mapping is best, do nothing.
    //But if an alternative mapping was best,
    //swap the coordinates according to the best mapping.
    if ( Eorig < minE ) {
      minE = Eorig;
    }
    else {
      if ( bestSwapID < 0 ) {
        cout << "Error, bestSwapID < 0" << endl;
        exit(0);
      }
      swapIndices( indices, altmaps[bestSwapID], swappedIndices );
      extractCoords( state1coords, swappedIndices, swappedCoords );
      setCoords( state1coords, indices, swappedCoords );
      cout << "--------" << endl;
      for ( size_t i = 0; i < indices.size(); i++ ) {
        cout << indices[i] << " -> " << swappedIndices[i] << endl;
      }
    }
}
