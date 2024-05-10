#include "PairTypeCutoffs.h"
#include <cstdlib>
#include <iostream>
#include <algorithm>

using namespace std;

PairTypeCutoffs::PairTypeCutoffs() :
  nAtomTypes(0),
  maxCutoff(0)
{
}

PairTypeCutoffs::~PairTypeCutoffs()
{
}

void PairTypeCutoffs::mapAtomsToTypes( const vector<int>& atomIndexToType_ ) {
  atomIndexToType = atomIndexToType_;
}

void PairTypeCutoffs::setTypeRadii( const vector<double>& radii ) {
  nAtomTypes = radii.size();
  //given N types,
  //setup the NxN lookup of cutoffs
  cutoffLookupStorage.resize( nAtomTypes*nAtomTypes );
  cutoffLookup.resize( nAtomTypes );
  for ( int i = 0; i < nAtomTypes; i++ ) {
    cutoffLookup[i] = &cutoffLookupStorage[i*nAtomTypes];
  }
  for ( int i = 0; i < nAtomTypes; i++ ) {
    for ( int j = i; j < nAtomTypes; j++ ) {
      cutoffLookup[i][j] = cutoffLookup[j][i] = radii[i] + radii[j];
    }
  }

  maxCutoff = *max_element( radii.begin(), radii.end() ) * 2.0;

}

void PairTypeCutoffs::setPairTypeCutoff( int atomtype1, int atomtype2, double cutoff ) {
  if ( atomtype1 >= nAtomTypes || atomtype2 >= nAtomTypes ) {
    cout << "Error in PairTypeCutoffs: atomtype exceeds nAtomTypes" << endl;
    exit(0);
  }
  cutoffLookup[atomtype1][atomtype2] =
    cutoffLookup[atomtype2][atomtype1] = cutoff;
  if ( cutoff > maxCutoff ) maxCutoff = cutoff;
}

void PairTypeCutoffs::check() {
  const size_t n = atomIndexToType.size();
  for ( size_t atom = 0; atom < n; atom++ ) {
    int atomtype = atomIndexToType[atom];
    if ( atomtype >= nAtomTypes ) {
      cout << "Error in PairTypeRepulsion : atomtype " << atomtype << " >= nAtomTypes " << nAtomTypes << endl;
      exit(0);
    }
  }
}
