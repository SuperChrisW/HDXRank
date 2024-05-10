#include "FrodaCutoffsBuilder.h"
#include "PairTypeCutoffs.h"
#include "ProteinInfo.h"
#include "NeighborTable.h"
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

FrodaCutoffsBuilder::FrodaCutoffsBuilder(
    const ProteinInfo &proteinInfo,
    const NeighborTable &nt,
    double mismatchTol_,
    double scaleFactor ) :
  polarGeometry(true),
  polarHRadius(0.5),
  vdwTol(scaleFactor),//0.85
  mismatchTol(mismatchTol_)
{
  pairTypeCutoffs = new PairTypeCutoffs;
  setupTypes();
  assignTypes( proteinInfo, nt );
  setupInteractionCutoffLookupTable();
  pairTypeCutoffs->check();
}

FrodaCutoffsBuilder::~FrodaCutoffsBuilder()
{
}

void FrodaCutoffsBuilder::setupTypes() {
  nameOfType.push_back( "H" );
  radiusOfType.push_back( 1.00 );
  chargeOfType.push_back( 0 );
  isPotentialDonorType.push_back( false );

  //The Hpolar type has a default radius of 1.00.
  //In certain situations it takes on the polarHRadius.
  nameOfType.push_back( "Hpolar" );
  radiusOfType.push_back( 1.00 );
  chargeOfType.push_back( 1 );
  isPotentialDonorType.push_back( false );

  nameOfType.push_back( "C" );
  radiusOfType.push_back( 1.70 );
  chargeOfType.push_back( 0 );
  isPotentialDonorType.push_back( false );

  nameOfType.push_back( "Cpolar" );
  radiusOfType.push_back( 1.70 );
  chargeOfType.push_back( 1 );
  isPotentialDonorType.push_back( false );

  nameOfType.push_back( "N" );
  radiusOfType.push_back( 1.55 );
  chargeOfType.push_back( -1 );
  isPotentialDonorType.push_back( false );

  nameOfType.push_back( "N_withAttachedH" );
  radiusOfType.push_back( 1.55 );
  chargeOfType.push_back( -1 );
  isPotentialDonorType.push_back( true );

  nameOfType.push_back( "O" );
  radiusOfType.push_back( 1.40 );
  chargeOfType.push_back( -1 );
  isPotentialDonorType.push_back( false );

  nameOfType.push_back( "O_withAttachedH" );
  radiusOfType.push_back( 1.40 );
  chargeOfType.push_back( -1 );
  isPotentialDonorType.push_back( true );

  nameOfType.push_back( "S" );
  radiusOfType.push_back( 1.80 );
  chargeOfType.push_back( 0 );
  isPotentialDonorType.push_back( false );

  nameOfType.push_back( "P" );
  radiusOfType.push_back( 1.80 );
  chargeOfType.push_back( 1 );
  isPotentialDonorType.push_back( false );

  nameOfType.push_back( "MG" );
  radiusOfType.push_back( 1.50 );
  chargeOfType.push_back( 0 );
  isPotentialDonorType.push_back( false );

  nameOfType.push_back( "MN" );
  radiusOfType.push_back( 1.40 );
  chargeOfType.push_back( 0 );
  isPotentialDonorType.push_back( false );

  nameOfType.push_back( "SI" );
  radiusOfType.push_back( 2.10 );
  chargeOfType.push_back( 1 );
  isPotentialDonorType.push_back( false );

  nameOfType.push_back( "FE" );
  radiusOfType.push_back( 1.50 );
  chargeOfType.push_back( 1 );
  isPotentialDonorType.push_back( false );

  nameOfType.push_back( "other" );
  radiusOfType.push_back( 1.50 );
  chargeOfType.push_back( 0 );
  isPotentialDonorType.push_back( false );

  nTypes = nameOfType.size();
  for ( int i=0; i<nTypes; i++ ) {
    mapNameToType[nameOfType[i]] = i;
  }

}

void FrodaCutoffsBuilder::assignTypes( const ProteinInfo &prot, const NeighborTable& nt ) {
  int n = prot.natoms();
  vector<int> atomtypes( n );

  int atomType;
  isHydrogen.resize(n, false);
  for ( int p=0; p<n; p++) {
    if ( prot.atom(p).elem() == "H" ) {
      isHydrogen[p] = true;
      // if any of p's neighbors is a N or O, atom is type Hpolar.
      // otherwise, it is type H
      atomType = mapNameToType["H"];
      vector<int>::const_iterator neighbor;
      for ( neighbor = nt[p].begin();
            neighbor != nt[p].end();
            neighbor++ ) {
        if ( prot.atom(*neighbor).elem() == "N" ||
             prot.atom(*neighbor).elem() == "O"){
          atomType = mapNameToType["Hpolar"];
          break;//break out of for loop.  Type has been set to Hpolar.
        }
      }
    }
    else if ( prot.atom(p).elem() == "C" ) {
      if ( prot.atom(p).name() == "C" ) atomType = mapNameToType["Cpolar"];
      else {
        int n_charged_neighbors = 0;
        // count the number of charged neighbors.
        vector<int>::const_iterator neighbor;
        for ( neighbor = nt[p].begin();
              neighbor != nt[p].end();
              neighbor++ ) {
          if ( prot.atom(*neighbor).elem() == "N" ||
               prot.atom(*neighbor).elem() == "O"){
              n_charged_neighbors++;
          }
        }
        if (n_charged_neighbors > 1){ // siteNumber carbon bonded to two Os or Ns gets siteNumber positive charge to balance them
          atomType = mapNameToType["Cpolar"];
        }
        else atomType = mapNameToType["C"];
      }
    }
    else if ( prot.atom(p).elem() == "N" ) {
      atomType = mapNameToType["N"];
      vector<int>::const_iterator neighbor;
      for ( neighbor = nt[p].begin();
            neighbor != nt[p].end();
            neighbor++ ) {
        if ( prot.atom(*neighbor).elem() == "H" ) {
          atomType = mapNameToType["N_withAttachedH"];
          break;
        }
      }
    }
    else if ( prot.atom(p).elem() == "O" ) {
      atomType = mapNameToType["O"];
      vector<int>::const_iterator neighbor;
      for ( neighbor = nt[p].begin();
            neighbor != nt[p].end();
            neighbor++ ) {
        if ( prot.atom(*neighbor).elem() == "H" ) {
          atomType = mapNameToType["O_withAttachedH"];
          break;
        }
      }
    }
    else if ( prot.atom(p).elem() == "S" ) atomType = mapNameToType["S"];
    else if ( prot.atom(p).elem() == "P" ) atomType = mapNameToType["P"];
    else if ( prot.atom(p).elem() == "MG" ) atomType = mapNameToType["MG"];
    else if ( prot.atom(p).elem() == "MN" ) atomType = mapNameToType["MN"];
    else if ( prot.atom(p).elem() == "SI" ) atomType = mapNameToType["SI"];
    else if ( prot.atom(p).elem() == "FE" ) atomType = mapNameToType["FE"];
    else {
      atomType = mapNameToType["other"];
      cout << "old-Froda-style Repulsion warning: atom " << p <<
              " has unrecognized element \"" << prot.atom(p).elem() << "\"" << endl;
    }

    atomtypes[p] = atomType;
  }

  pairTypeCutoffs->mapAtomsToTypes( atomtypes );

}

void FrodaCutoffsBuilder::setupInteractionCutoffLookupTable() {
  //initialize cutoffs according to radii
  pairTypeCutoffs->setTypeRadii( radiusOfType );
  ofstream outfile("frodacutoffs.txt", ios::out );

  //overwrite cutoffs according to the original FRODA logic
  for ( int t1 = 0; t1 < nTypes; t1++ ) {
    for ( int t2 = t1; t2 < nTypes; t2++ ) {
      double cutoff = calcInteractionCutoffForTypePair( t1, t2 ) + mismatchTol;
      pairTypeCutoffs->setPairTypeCutoff( t1, t2, cutoff );
      outfile << t1 << " " << t2 << " " << cutoff << '\n';
    }
  }
  outfile.close();
}

double FrodaCutoffsBuilder::calcInteractionCutoffForTypePair( int t1, int t2 ) const {
  double r1 = radiusOfType[t1];
  double r2 = radiusOfType[t2];
  double cutoff;
  if ( polarGeometry && (chargeOfType[t1] * chargeOfType[t2]) < 0  ) {
    if ( nameOfType[t1] == "Hpolar" ) cutoff = (polarHRadius + r2) * vdwTol;
    else if ( nameOfType[t2] == "Hpolar" ) cutoff = (r1 + polarHRadius) * vdwTol;
    else cutoff = ( r1 + r2 ) * vdwTol * 0.92; //tuned to backbone C O in alpha helix
  }
  else if ( polarGeometry &&
        (( isPotentialDonorType[t1] && chargeOfType[t2] == -1 ) ||
         ( isPotentialDonorType[t2] && chargeOfType[t1] == -1 )) ) {
    // we assume that the pair is a donor-acceptor pair
    // participating in a hydrogen bond.
    cutoff = ( r1 + r2 ) * vdwTol * 0.92;
  }
  else cutoff = ( r1 + r2 ) * vdwTol;

  return cutoff;

}
