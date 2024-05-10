/*
 * SideChainTorsion.cpp
 *
 *  Created on: Nov 14, 2009
 *      Author: dwfarrel
 */

#include "SideChainTorsion.h"
#include "ProteinInfo.h"
#include "NeighborTable.h"
#include "RigidUnitSystem.h"
#include <cmath>
#include <string>

using namespace std;

SideChainTorsionInitializer::SideChainTorsionInitializer(
    const ProteinInfo& prot_,
    const NeighborTable& nt_,
    const std::vector<Vec3>& coords_,
    SideChainTorsionContainer& tor_ ) :
  prot( prot_ ),
  nt( nt_ ),
  coords( coords_ ),
  tor( tor_ ),
  k( 10.0 ) {

  //only consider atoms within the standard amino acids,
  //excluding proline and glycine.  The 18 remaining amino
  //acids all have at least one rotatable bond subject to
  //a torsion constraint.
  resList.insert( "LEU" );
  resList.insert( "ALA" );
  resList.insert( "VAL" );
  resList.insert( "MET" );
  resList.insert( "ILE" );
  resList.insert( "TRP" );
  resList.insert( "PHE" );
  resList.insert( "SER" );
  resList.insert( "THR" );
  resList.insert( "CYS" );
  resList.insert( "ASN" );
  resList.insert( "GLN" );
  resList.insert( "TYR" );
  resList.insert( "HIS" );
  resList.insert( "ASP" );
  resList.insert( "GLU" );
  resList.insert( "LYS" );
  resList.insert( "ARG" );
}

SideChainTorsionInitializer::~SideChainTorsionInitializer() {
}

void SideChainTorsionInitializer::setupConstraints() {

  double thetadeg = 55;
  double thetarad = thetadeg*M_PI/180.0;
  double costheta = cos(thetarad);
  double sintheta = sin(thetarad);

  int natoms = prot.natoms();
  const ProteinInfo::Atoms& atoms( prot.atoms() );

  vector<char> is_CNS_SP3( natoms, 0 );
  for ( int i = 0; i < natoms; i++ ) {
    //only consider atoms from the allowed list of residues
    if ( resList.find( atoms[i].resi().name() ) == resList.end() ) continue;
    if ( nt[i].size() == 4 &&
         ( atoms[i].elem() == "C" || atoms[i].elem() == "N" || atoms[i].elem() == "S" ) ) {
      is_CNS_SP3[i] = 1;
    }
  }

  //NOTE:  To identify torsion constraints, the loop below looks for
  //a sequence of bonded atoms (m -> i -> j -> n) such that i and j
  //each have 4 neighbors, and i and j are either C, N, or S.
  //This logic could fail if given a tetrahedrally-connected ring like proline, because
  //it would create constraints within the ring which is rigid.  The logic could not fail
  //for flat rings like phenol, because the atoms do not have 4 neighbors.
  //It could also potentially fail for any arrangement in which atoms
  //m and n are in the same rigid unit.
  //Because of this, I am using a specific residue list to exclude proline.
  for ( int i = 0; i < natoms; i++ ) {
    if ( !is_CNS_SP3[i] ) continue;
    const vector<int>* neighbors;
    vector<int>::const_iterator it;
    neighbors = &nt[i];
    for ( it = neighbors->begin(); it != neighbors->end(); it++ ) {
      int j = *it;
      if ( i >= j ) continue;
      if ( !is_CNS_SP3[j] ) continue;

      //if we get here, then both i and j are bonded and are both SP3 C, N, or S
      //So, we want to make 1-4 cutoffs between the neighbors of i and the
      //neighbors of j.

      //go over each neighbor of i, and each neighbor of j,
      //excluding of course i and j themselves.
      const vector<int>* neighbors1;
      vector<int>::const_iterator it1;
      const vector<int>* neighbors2;
      vector<int>::const_iterator it2;
      neighbors1 = &nt[i];
      neighbors2 = &nt[j];
      int m;
      int n;
      for ( it1 = neighbors1->begin(); it1 != neighbors1->end(); it1++ ) {
        m = *it1;
        if ( m == j ) continue;
        for ( it2 = neighbors2->begin(); it2 != neighbors2->end(); it2++ ) {
          n = *it2;
          if ( n == i ) continue;
          //we have a 1-4 relationship between m and n, across two sp3 atoms i and j
          if ( m == n ||
               nt.isFirstNeighbor( m, n ) ||
               nt.isSecondNeighbor( m, n ) ) continue;

          double constraintdist = calc14Constraint( m, i, j, n, costheta, sintheta );
          int p1 = min( m, n );
          int p4 = max( m, n );
          tor.insert( p1, p4, MinDistConstraint( &coords, k, p1, p4, constraintdist ) );
        }
      }
    }
  }
}

void SideChainTorsionInitializer::excludeRigidPairs( const RigidUnitSystem* sys ) {
  SideChainTorsionContainer::iterator it = tor.begin();
  while ( it != tor.end() ) {
    if ( sys->doPointsBelongToSameRigidUnit( it->getp1(), it->getp2() ) ) {
      it = tor.erase( it );
    }
    else it++;
  }
}

double SideChainTorsionInitializer::calc14Constraint(
    int a1, int a2, int a3, int a4,
    double costheta, double sintheta ) {

  Vec3 A42 = coords[a4] - coords[a2];
  Vec3 A32 = coords[a3] - coords[a2];
  Vec3 A12 = coords[a1] - coords[a2];
  double A32_length = sqrt( A32.norm2() );
  Vec3 axis = A32 * (1/A32_length);

  Vec3 A1;
  A1.z = A12.dot( axis );
  A1.x = sqrt( A12.norm2() - A1.z*A1.z );
  A1.y = 0;

  Vec3 A4;
  A4.z = A42.dot( axis );
  double A4_perp = sqrt( A42.norm2() - A4.z*A4.z );
  A4.x = A4_perp*costheta;
  A4.y = A4_perp*sintheta;

  //Now calculate the 1-4 distance of this arrangement.
  return sqrt( A1.dist2( A4 ) );
}
