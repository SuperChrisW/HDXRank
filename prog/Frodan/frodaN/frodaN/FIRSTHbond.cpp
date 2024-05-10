/*
 * FIRSTHbond.cpp
 *
 *  Created on: Mar 10, 2009
 *      Author: dwfarrel
 */

#include "FIRSTHbond.h"
#include "ProteinInfo.h"
#include <cmath>

using namespace std;

FIRSTHbond::FIRSTHbond(
    const ProteinInfo& proteinInfo_,
    const NeighborTable& nt_,
    const std::vector<Vec3>& coords_ ) :
  proteinInfo( proteinInfo_ ),
  nt( nt_ ),
  coords( coords_ ),
  lookupHBAtomInfo( proteinInfo_.natoms() )
{
  SB_hyd_accpt_maxdist = 4.0;
  SB_donor_accpt_maxdist = 5.0;
  SB_donor_hyd_accpt_minangle = 100.0;
  SB_hyd_accpt_base_minangle = 80.0;

  HB_hyd_accpt_maxdist = 4.0;
  HB_donor_accpt_maxdist = 5.0;
  HB_donor_hyd_accpt_minangle = 100.0;

  cutoff_SR_normal_plane_angle = 30.0;
  cutoff_SR_center_normal_angle = 40.0;
  cutoff_SR_distance = 5.5;
  CNmaxDoubleBondDistance = 1.4;

  SB_hyd_accpt_maxdist2 = SB_hyd_accpt_maxdist*SB_hyd_accpt_maxdist;
  SB_donor_accpt_maxdist2 = SB_donor_accpt_maxdist*SB_donor_accpt_maxdist;
  SB_donor_hyd_accpt_maxcosangle = cos( SB_donor_hyd_accpt_minangle*M_PI/180.0 );
  SB_hyd_accpt_base_maxcosangle = cos( SB_hyd_accpt_base_minangle*M_PI/180.0 );

  HB_hyd_accpt_maxdist2 = HB_hyd_accpt_maxdist*HB_hyd_accpt_maxdist;
  HB_donor_accpt_maxdist2 = HB_donor_accpt_maxdist*HB_donor_accpt_maxdist;
  HB_donor_hyd_accpt_maxcosangle = cos( HB_donor_hyd_accpt_minangle*M_PI/180.0 );

  CNmaxDoubleBondDistance2 = CNmaxDoubleBondDistance*CNmaxDoubleBondDistance;

  assignHydrogenBondStatus();
}

FIRSTHbond::~FIRSTHbond() {
  // TODO Auto-generated destructor stub
}

void FIRSTHbond::assignHydrogenBondStatus(){
  int N = proteinInfo.natoms();
  for ( int i = 0; i < N; i++ ) {
    if ( nt[i].size() == 0 ) continue;
    if ( proteinInfo.atom(i).elem() == "N" ) {
      lookupHBAtomInfo[i] = determineNitrogenHydrogenBondStatus( i );
    }
    else if( proteinInfo.atom(i).elem() == "O" ) {
      lookupHBAtomInfo[i] = determineOxygenHydrogenBondStatus( i );
    }
    else if( proteinInfo.atom(i).elem() == "S" ) {
      lookupHBAtomInfo[i] = determineSulfurHydrogenBondStatus( i );
    }
  }
  for ( int i = 0; i < N; i++ ) {
    if ( proteinInfo.atom(i).elem() == "H" &&
         nt[i].size() == 1 &&
         lookupHBAtomInfo[ nt[i][0] ].isDonor )
      lookupHBAtomInfo[i].isDonorH = 1;
  }
}

HBAtomInfo FIRSTHbond::determineNitrogenHydrogenBondStatus( int i ) {
  HBAtomInfo info;

  //regular backbone Nitrogens
  if ( proteinInfo.atom(i).name() == "N" && nt[i].size() == 3 ) {
    info.isDonor = 1;
    info.isSP2 = 1;
    return info;
  }

  //Terminal backbone Nitrogens, and any other Nitrogen with 4 neighbors
  if ( nt[i].size() == 4 ) {
    info.isDonor = 1;
    info.isSP3 = 1;
    info.isCharged = 1;
    return info;
  }

  if ( proteinInfo.atom(i).resi().name() == "ARG" ) {
    info.isDonor = 1;
    info.isSP2 = 1;
    info.isCharged = 1;
    return info;
  }

  if ( proteinInfo.atom(i).resi().name() == "LYS" ) {
    info.isDonor = 1;
    info.isSP3 = 1;
    info.isCharged = 1;
    return info;
  }

  if( proteinInfo.atom(i).resi().name() == "ASN" || proteinInfo.atom(i).resi().name() == "GLN" ) {
    info.isDonor = 1;
    info.isSP2 = 1;
    return info;
  }

  if( proteinInfo.atom(i).resi().name() == "HIS" || proteinInfo.atom(i).resi().name() == "HID" ||
      proteinInfo.atom(i).resi().name() == "HIE" ) {

    info.isSP2 = 1;
    info.isCharged = 1;

    //determine whether donor or acceptor
    bool hasHneighbor = false;
    int Nneigh = nt[i].size();
    for ( int neigh = 0; neigh < Nneigh; neigh++ ) {
      int j = nt[i][neigh];
      if ( proteinInfo.atom(j).elem() == "H" ) {
        hasHneighbor = true;
        break;
      }
    }
    if ( hasHneighbor ) info.isDonor = 1;
    else info.isAcceptor = 1;
    return info;
  }

  if( proteinInfo.atom(i).resi().name() == "TRP" ) {
    info.isDonor = 1;
    info.isSP2 = 1;
    return info;
  }

  // If the nitrogen does not belong to a standard amino acid, set the
  // hydrogen bond status based on the total number of bonded neighbors
  // and the total number of neighbors that are nitrogen.
  // BMH 10.26.04 These default rules check only the valency to assign
  //   hybridization. Really should check the bond order of each bond
  //   to accurately determine the hybridization.
  //////////////////////////////////////////////////////////////////////
  if( nt[i].size() == 2 ) {
    info.isAcceptor = 1;
    info.isSP2 = 1;

    //determine charged status
    int neigh0 = nt[i][0];
    int neigh1 = nt[i][1];
    if ( !( proteinInfo.atom(neigh0).elem() == "C" && coords[i].dist2( coords[neigh0] ) <= CNmaxDoubleBondDistance2 ||
            proteinInfo.atom(neigh1).elem() == "C" && coords[i].dist2( coords[neigh1] ) <= CNmaxDoubleBondDistance2 ) )
    {
      info.isCharged = 1;
    }
    return info;
  }

  if( nt[i].size() == 3 ) {
    info.isDonor = 1;
    info.isSP2 = 1;
    return info;
  }

  if ( nt[i].size() >= 5 ) {
    cout << "Nitrogen " << i << " has 5 or more covalently bonded neighbors." << endl;
    exit(0);
  }

  // Default. Return Donor SP2.
  //////////////////////////////////////////////////////////////////////

  info.isDonor = 1;
  info.isSP2 = 1;
  return info;
}

HBAtomInfo FIRSTHbond::determineOxygenHydrogenBondStatus( int i ) {
  HBAtomInfo info;
  //backbone
  if( proteinInfo.atom(i).name() == "O" ) {
    info.isAcceptor = 1;
    info.isSP2 = 1;
    return info;
  }

  //backbone terminal carboxyl
  if ( nt[i].size() == 1 &&
        proteinInfo.atom( nt[i][0] ).name() == "C" ) {
    info.isAcceptor = 1;
    info.isSP2 = 1;
    info.isCharged = 1;
    return info;
  }

  if ( proteinInfo.atom(i).resi().name() == "ASP" || proteinInfo.atom(i).resi().name() == "GLU" ) {
    info.isAcceptor = 1;
    info.isSP2 = 1;
    info.isCharged = 1;
    return info;
  }

  if( proteinInfo.atom(i).resi().name() == "SER" || proteinInfo.atom(i).resi().name() == "THR" ) {
    info.isDonor = 1;
    info.isAcceptor = 1;
    info.isSP3 = 1;
    return info;
  }

  if( proteinInfo.atom(i).resi().name() == "TYR" ) {
    info.isDonor = 1;
    info.isAcceptor = 1;
    info.isSP2 = 1;
    return info;
  }

  if( proteinInfo.atom(i).resi().name() == "ASN" || proteinInfo.atom(i).resi().name() == "GLN" ) {
    info.isAcceptor = 1;
    info.isSP2 = 1;
    return info;
  }

  // if the oxygen has only one covalent bond, see if it's part of a
  // carboxylate, or not.
  //////////////////////////////////////////////////////////////////////
  if( nt[i].size() == 1 ) { // If the current oxygen has one neighbor, check
    info.isAcceptor = 1;
    info.isSP2 = 1;

    //determine charged status
    int carbon = nt[i][0];
    int numOxygens = 0;
    int Nneigh = nt[carbon].size();
    for ( int neigh = 0; neigh < Nneigh; neigh++ ) {
      int j = nt[carbon][neigh];
      if ( proteinInfo.atom(j).elem() == "O" ) {
        numOxygens++;
      }
    }
    if ( numOxygens == 2 ) {
      info.isCharged = 1;
    }
    return info;
  }

  // if the oxygen has two neighbors, one of which is a hydrogen, set it
  // to donor and acceptor.
  //////////////////////////////////////////////////////////////////////
  if( nt[i].size() == 2 ){
    info.isAcceptor = 1;
    info.isSP3 = 1;

    if( proteinInfo.atom( nt[i][0] ).elem() == "H" ||
        proteinInfo.atom( nt[i][1] ).elem() == "H" ) {
      info.isDonor = 1;
    }
    return info;
  }

  // If the oxygen is isolated.
  //////////////////////////////////////////////////////////////////////
  if( nt[i].size() == 0 ){
    return info;
  }

  // If the oxygen has 3 covalent bonds, something went wrong.
  //////////////////////////////////////////////////////////////////////
  cout << "WARNING: Oxygen " << i << " has " << nt[i].size() << " covalent neighbors" << endl;
  return info;
}

HBAtomInfo FIRSTHbond::determineSulfurHydrogenBondStatus( int i ) {
  HBAtomInfo info;
  if ( proteinInfo.atom(i).resi().name() == "MET" ) {
    info.isAcceptor = 1;
    info.isSP3 = 1;
    return info;
  }

  // If it's a Cysteine residue, mark as an sp3 donor-acceptor. If it's
  // a cystine involved in a disulfide, mark it as an sp3 acceptor only.
  //////////////////////////////////////////////////////////////////////
  if ( proteinInfo.atom(i).resi().name() == "CYS" ||
       proteinInfo.atom(i).resi().name() == "CYX" ) {
    info.isAcceptor = 1;
    info.isSP3 = 1;

    //determine donor status
    bool hasHneighbor = false;
    int Nneigh = nt[i].size();
    for ( int neigh = 0; neigh < Nneigh; neigh++ ) {
      int j = nt[i][neigh];
      if ( proteinInfo.atom(j).elem() == "H" ) {
        hasHneighbor = true;
        break;
      }
    }
    if ( hasHneighbor ) {
      info.isDonor = 1;
    }
    return info;
  }

  // Process sulfur atom in an unknown group.
  //////////////////////////////////////////////////////////////////////
  if ( nt[i].size() == 1 ) {
    info.isAcceptor = 1;
    info.isSP2 = 1;
    info.isCharged = 1;
    return info;
  }
  if ( nt[i].size() == 2 ) {
    info.isAcceptor = 1;
    info.isSP3 = 1;

    if ( proteinInfo.atom( nt[i][0] ).elem() == "H" ||
         proteinInfo.atom( nt[i][1] ).elem() == "H" ) {
      info.isDonor = 1;
    }
    return info;
  }
  if ( nt[i].size() == 3 ) {
    info.isDonor = 1;
    info.isSP3 = 1;
    info.isCharged = 1;
    return info;
  }

  cout << "Warning: Could not determine hydrogen bond status for sulfur atom " << i << endl;
  return info;
}

double FIRSTHbond::energy( int h, int a ) const {
  if ( !lookupHBAtomInfo[h].isDonorH || !lookupHBAtomInfo[a].isAcceptor ) return 0;
  if ( nt.isThirdNeighbor( h, a ) ||
       nt.isSecondNeighbor( h, a ) ||
       nt.isFirstNeighbor( h, a ) ) return 0;

  double e = 0;
  //the fact that isDonorH is true means that the h does have 1 neighbor,
  //and that neighbor is a donor.
  int d = nt[h][0];
  if( lookupHBAtomInfo[d].isCharged && lookupHBAtomInfo[a].isCharged )
    e = saltBridgeEnergy(d, h, a);
  else e = hydrogenBondEnergy(d, h, a);
  return ( e > 0 ) ? 0 : e;
}

bool FIRSTHbond::isDonor( int i ) const {
  return lookupHBAtomInfo[i].isDonor;
}

bool FIRSTHbond::isAcceptor( int i ) const {
  return lookupHBAtomInfo[i].isAcceptor;
}

bool FIRSTHbond::isDonorH( int i ) const {
  return lookupHBAtomInfo[i].isDonorH;
}

float FIRSTHbond::hydrogenBondEnergy( int don, int hyd, int acc ) const {
  // Make sure the pair satisfy the geometric cutoff.
  float DA_distance2 = coords[don].dist2( coords[acc] );
  if ( DA_distance2 > HB_donor_accpt_maxdist2 )
    return 0;
  if ( coords[hyd].dist2( coords[acc] ) > HB_hyd_accpt_maxdist2 )
    return 0;
  if ( computeCosAngle(don, hyd, acc) > HB_donor_hyd_accpt_maxcosangle )
    return 0;

  // Compute the distance dependent term of the total energy.
  float R_s = 2.8;

  float ratio =  R_s / sqrt(DA_distance2);
  float ratio2 = ratio*ratio;
  float ratio3 = ratio2*ratio;
  float ratio5 = ratio3*ratio2;
  float ratio6 = ratio3*ratio3;
  float E_distance = 8.0 * ( 5*ratio6*ratio6 - 6*ratio5*ratio5 );

  // 4. Compute the angular prefactor [cos^2(theta) * e^(-(pi - theta)^6)].
  //    All of the hybridization-dependent angular terms contain this prefactor.
  //////////////////////////////////////////////////////////////////////
  float costheta = computeCosAngle( don, hyd, acc );
  float piMinusTheta = M_PI - acos(costheta);
  float piMinusTheta2 = piMinusTheta*piMinusTheta;
  float piMinusTheta3 = piMinusTheta2*piMinusTheta;
  float prefactor = costheta*costheta*exp( -piMinusTheta3*piMinusTheta3 );

  // 6. Compute the additional angular terms, as necessary, for the 4 donor-
  //    acceptor hybridization pairs.
  ////////////////////////////////////////////////////////////////////////////////

  float E_angular;

  // SP2 donor - SP2 acceptor hydrogen bond.
  if( lookupHBAtomInfo[don].isSP2 && lookupHBAtomInfo[acc].isSP2 ){

    vector<Vec3> hydNormals;
    findNormals( hyd, &hydNormals );

    vector<Vec3> accNormals;
    findNormals( acc, &accNormals );

    float lowestCosAngle = 1;
    int Nneigh = nt[acc].size();
    for( int neigh = 0; neigh < Nneigh; neigh++ ){
      int base = nt[acc][neigh];
      if ( proteinInfo.atom(base).elem() == "H" ) continue; //not sure why we're excluding this case
      float cosphi = computeCosAngle(hyd, acc, base);
      if ( cosphi < lowestCosAngle ) lowestCosAngle = cosphi;
      if( cosphi > 0 )
        return 0;
    }

    int nHydNormals = hydNormals.size();
    int nAccNormals = accNormals.size();
    for( int i = 0; i < nHydNormals; i++ ) {
      for( int j = 0; j < nAccNormals; j++ ) {
        float cosgamma = computeOutOfPlaneCosAngle( hydNormals[i], accNormals[j] );
        if ( cosgamma < lowestCosAngle ) lowestCosAngle = cosgamma;
      }
    }
    E_angular = lowestCosAngle*lowestCosAngle*prefactor;
  }

  // SP3 donor - SP2 acceptor hydrogen bond.
  //Note - replace with loop over all non-H neighbors of the O acceptor
  else if ( lookupHBAtomInfo[don].isSP3 && lookupHBAtomInfo[acc].isSP2 ) {
    float lowestCosPhi = 1;
    int Nneigh = nt[acc].size();
    for( int neigh = 0; neigh < Nneigh; neigh++ ){
      int base = nt[acc][neigh];
      if ( proteinInfo.atom(base).elem() == "H" ) continue; //not sure why we're excluding this case
      float cosphi = computeCosAngle(hyd, acc, base);
      if ( cosphi < lowestCosPhi ) lowestCosPhi = cosphi;
      if( cosphi > 0 )
        return 0;
    }

    E_angular = prefactor*lowestCosPhi*lowestCosPhi;
  }

  // SP2 donor - SP3 acceptor hydrogen bond.
  else if ( lookupHBAtomInfo[don].isSP2 && lookupHBAtomInfo[acc].isSP3 ) {
    E_angular = prefactor*prefactor;
  }

  // SP3 donor - SP3 acceptor hydrogen bond.
  //Note - replace with loop over all non-H neighbors of the O acceptor, or
  //       if water, it is both H's.
  else if ( lookupHBAtomInfo[don].isSP3 && lookupHBAtomInfo[acc].isSP3 ) {
    float lowestPhi = M_PI;
    int Nneigh = nt[acc].size();
    for( int neigh = 0; neigh < Nneigh; neigh++ ){
      int base = nt[acc][neigh];

      if ( proteinInfo.atom(base).elem() == "H" ) continue; //not sure why we're excluding this case

      float phi = acos( computeCosAngle(hyd, acc, base) );
      if ( phi < lowestPhi ) lowestPhi = phi;
    }
    const float phi0 = 109.50*M_PI/180.0;
    float cosAngleDiff = cos(lowestPhi - phi0);
    E_angular = prefactor*cosAngleDiff*cosAngleDiff;
  }
  else return 0;

  return E_distance*E_angular;
}

float FIRSTHbond::computeCosAngle( int i1, int i2, int i3 ) const {
  Vec3 r21 = coords[i1] - coords[i2];
  Vec3 r23 = coords[i3] - coords[i2];
  return r21.dot( r23 )/sqrt( r21.norm2()*r23.norm2() );
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   This function is used when computing the hydrogen bond energy when both
//   the donor and acceptor are sp2 hybridized. The neighbors of sp2 hybridized
//   atoms lie in a plane. The function calculates the normals from each sp2
//   plane, and returns the angle between the normals.
////////////////////////////////////////////////////////////////////////////////
float FIRSTHbond::computeOutOfPlaneCosAngle( const Vec3& unitnormal1, const Vec3& unitnormal2 ) const {
  float cosgamma ( unitnormal1.dot( unitnormal2 ) );
  return -abs( cosgamma );
}

void FIRSTHbond::getUnitNormalToPlane( int i1, int i2, const int i3, Vec3* c) const {

  Vec3 a = coords[i2] - coords[i1];
  Vec3 b = coords[i3] - coords[i1];

  c->x = a.y*b.z - a.z*b.y;
  c->y = a.z*b.x - a.x*b.z;
  c->z = a.x*b.y - a.y*b.x;

  *c /= sqrt( c->norm2() );
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Find a set of three consecutively bonded atoms centered on the input atom.
//   These three atoms will form an angle, which is used in the energy functions.
//   If the input atom is singly coordinated, search for all triples in which the
//   the input atom is a terminal vertex on the two-edge chain.
////////////////////////////////////////////////////////////////////////////////
void FIRSTHbond::findNormals( int site1, vector<Vec3>* normals  ) const {

  int site2;
  int site3;

  normals->clear();
  Vec3 normal;

  if ( nt[site1].size() == 1 ) {
    site2 = nt[site1][0];
    int Nneigh = nt[site2].size();
    for( int neigh = 0; neigh < Nneigh; neigh++ ){
      site3 = nt[site2][neigh];
      if ( site3 == site1 ) continue;
      getUnitNormalToPlane( site1, site2, site3, &normal );
      normals->push_back( normal );
    }
  }
  else {
    int Nneigh = nt[site1].size();
    for ( int neigh1 = 0; neigh1 < Nneigh; neigh1++ ) {
      site2 = nt[site1][neigh1];
      for ( int neigh2 = neigh1 + 1; neigh2 < Nneigh; neigh2++ ) {
        site3 = nt[site1][neigh2];
        getUnitNormalToPlane( site1, site2, site3, &normal );
        normals->push_back( normal );
      }
    }
  }

}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Compute the energy for salt bridges.
//
//   The salt bridge function has the following form: (L-J 10-12 type)
//
//                 {   ( R_s )12      ( R_s )10 }
//     E_sb = 10.0*{ 5*(-----)   -  6*(-----)   }
//                 (   (R + a)        (R + a)   }
//
//   where: R = donor-acceptor distance.
//          The prefactor 10.0 is the well depth.
//          R_s is the donor-acceptor distance at the minimum.
////////////////////////////////////////////////////////////////////////////////
float FIRSTHbond::saltBridgeEnergy( int don, int hyd, int acc ) const {
  // Make sure the pair satisfy the geometric cutoff.
  float DA_distance2 = coords[don].dist2( coords[acc] );

  if ( DA_distance2 > SB_donor_accpt_maxdist2 )
    return 0;
  if ( coords[hyd].dist2( coords[acc] ) > SB_hyd_accpt_maxdist2 )
    return 0;
  if( computeCosAngle(don, hyd, acc) > SB_donor_hyd_accpt_maxcosangle )
    return 0;

  int Nneigh = nt[acc].size();
  for( int neigh = 0; neigh < Nneigh; neigh++ ){
    int base = nt[acc][neigh];
    if( computeCosAngle(hyd, acc, base) > SB_hyd_accpt_base_maxcosangle )
      return 0;
  }

  float R_s = 3.2;
  float a = 0.375;

  float ratio = (R_s / ( sqrt(DA_distance2) + a) );
  float ratio2 = ratio*ratio;
  float ratio3 = ratio2*ratio;
  float ratio5 = ratio3*ratio2;
  float ratio6 = ratio3*ratio3;

  return 10.0 * ( 5*ratio6*ratio6 - 6*ratio5*ratio5 );
}
