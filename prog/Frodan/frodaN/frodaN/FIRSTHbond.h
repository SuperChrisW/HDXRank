/*
 * FIRSTHbond.h
 *
 *  Created on: Mar 10, 2009
 *      Author: dwfarrel
 */

#ifndef FIRSTHBOND_H_
#define FIRSTHBOND_H_

#include "NeighborTable.h"
#include "Vec3.h"
#include <vector>
#include <string>
class ProteinInfo;

class HBAtomInfo {
public:
  HBAtomInfo() : isDonor(0), isDonorH(0), isAcceptor(0), isCharged(0), isSP2(0), isSP3(0) {}
  ~HBAtomInfo() {}
  void clear() { isDonor = 0; isDonorH = 0; isAcceptor = 0; isCharged = 0; isSP2 = 0; isSP3 = 0; }
  char isDonor;
  char isDonorH;
  char isAcceptor;
  char isCharged;
  char isSP2;
  char isSP3;
};

class FIRSTHbond {
public:
  FIRSTHbond(
    const ProteinInfo& proteinInfo,
    const NeighborTable& nt,
    const std::vector<Vec3>& coords );

  virtual ~FIRSTHbond();

  bool isDonor( int i ) const;
  bool isDonorH( int i ) const;
  bool isAcceptor( int i ) const;
  double energy( int h, int a ) const;

private:
  const ProteinInfo& proteinInfo;
  const NeighborTable& nt;
  const std::vector<Vec3>& coords;
  std::vector< HBAtomInfo > lookupHBAtomInfo;

  float SB_hyd_accpt_maxdist;
  float SB_donor_accpt_maxdist;
  float SB_donor_hyd_accpt_minangle;
  float SB_hyd_accpt_base_minangle;

  float HB_hyd_accpt_maxdist;
  float HB_donor_accpt_maxdist;
  float HB_donor_hyd_accpt_minangle;

  float cutoff_SR_normal_plane_angle;
  float cutoff_SR_center_normal_angle;
  float cutoff_SR_distance;

  float CNmaxDoubleBondDistance;

  float SB_hyd_accpt_maxdist2;
  float SB_donor_accpt_maxdist2;
  float SB_donor_hyd_accpt_maxcosangle;
  float SB_hyd_accpt_base_maxcosangle;

  float HB_hyd_accpt_maxdist2;
  float HB_donor_accpt_maxdist2;
  float HB_donor_hyd_accpt_maxcosangle;

  float CNmaxDoubleBondDistance2;
  void assignHydrogenBondStatus();
  HBAtomInfo determineNitrogenHydrogenBondStatus( int i );
  HBAtomInfo determineOxygenHydrogenBondStatus( int i );
  HBAtomInfo determineSulfurHydrogenBondStatus( int i );
  float hydrogenBondEnergy( int don, int hyd, int acc ) const;
  float computeCosAngle( int i1, int i2, int i3 ) const;
  float computeOutOfPlaneCosAngle( const Vec3& unitnormal1, const Vec3& unitnormal2 ) const;
  void getUnitNormalToPlane( int i1, int i2, const int i3, Vec3* c ) const;
  void findNormals( int site1, std::vector<Vec3>* normals ) const;
  float saltBridgeEnergy( int don, int hyd, int acc ) const;

};

#endif /* FIRSTHBOND_H_ */
