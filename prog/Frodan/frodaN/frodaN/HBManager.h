/*
 * HBManager.h
 *
 *  Created on: Jun 11, 2009
 *      Author: dwfarrel
 */

#ifndef HBMANAGER_H_
#define HBMANAGER_H_

class FIRSTHbond;
#include "VL.h"
#include "HBConstraint.h"
#include "NeighborTable.h"
#include <vector>

class HAPair {
public:
  HAPair() {}
  HAPair( int h_, int a_ ) : h(h_), a(a_) {}
  ~HAPair() {}
  int h;
  int a;
};

class VerletCollectorHAPairs {
public:
  VerletCollectorHAPairs(
    const std::vector<char>& isDonorH_, const std::vector<char>& isAcceptor_ ) :
      isDonorH( isDonorH_ ), isAcceptor( isAcceptor_) {}
  ~VerletCollectorHAPairs() {}
  void clearPairs() { pairs.clear(); }
  void addPair( int p1, int p2, double cutoff ) {
    if ( isDonorH[p1] && isAcceptor[p2] )
      pairs.push_back( HAPair(p1,p2) );
    else if ( isDonorH[p2] && isAcceptor[p1] )
      pairs.push_back( HAPair(p2,p1) );
  }
  typedef std::vector< HAPair >::const_iterator const_iterator;
  const_iterator begin() const { return pairs.begin(); }
  const_iterator end() const { return pairs.end(); }
private:
  std::vector< HAPair > pairs;
  const std::vector<char>& isDonorH;
  const std::vector<char>& isAcceptor;
};

class VerletCutoffsHAPairs {
public:
  VerletCutoffsHAPairs( const std::vector<char>& isDonorH_, const std::vector<char>& isAcceptor_, const NeighborTable* nt_ ) :
    isDonorH( isDonorH_ ), isAcceptor( isAcceptor_ ), nt( nt_ ) {}
  ~VerletCutoffsHAPairs() {}
  void getCutoff( int p1, int p2, bool& isPairExcluded, double& cutoff ) const {
    cutoff = 5.0;
    isPairExcluded =
      !( ( isDonorH[p1] && isAcceptor[p2] ) || ( isDonorH[p2] && isAcceptor[p1] ) ) ||
      nt->isFirstNeighbor( p1, p2) ||
      nt->isSecondNeighbor( p1, p2 ) ||
      nt->isThirdNeighbor( p1, p2 );// ||
      //proteinInfo->inSameResidue( p1, p2 );
  }
  double getMaxCutoff() const { return 5.0; }
private:
  const std::vector<char>& isDonorH;
  const std::vector<char>& isAcceptor;
  const NeighborTable* nt;
};

class ProteinInfo;
#include "Vec3.h"
class ForbidList;

//Note, all variables that are input in the constructor
//must remain valid throughout the existence of the HBManager object.
class HBManager
{
public:
  HBManager(
    const ProteinInfo& proteinInfo,
    const NeighborTable& nt,
    const std::vector<Vec3>& coords,
    BBHBContainer& bb_hbonds,
    HBContainer& hbonds );
  virtual ~HBManager();

  void disableAnglesInNewConstraints() { doAngleConstraints = false; }
  void enableAnglesInNewConstraints() { doAngleConstraints = true; }
  void enableCarefulChecking() { carefulChecking = true; }
  void disableCarefulChecking() { carefulChecking = false; }
  void setEcutoff( double E ) { Ecutoff = E; }

  bool isPairExcluded( int p1, int p2 ) const;
  void addConstraint( int h, int a, bool& success, bool breakable = true );
  void tighten();
  void findnew( bool breakable = true );

  void resetCounts() { nAdded_ = nRemoved_ = 0; }
  int nAdded() const { return nAdded_; }
  int nRemoved() const { return nRemoved_; }
  size_t size() const { return bb_hbonds.size() + hbonds.size(); }

  void breakAllBreakable();
  void attachForbidList( ForbidList* forbidlist_ ) { forbidlist = forbidlist_; }

private:
  const NeighborTable& nt;
  const std::vector<Vec3>& coords;
  BBHBContainer& bb_hbonds;
  HBContainer& hbonds;
  FIRSTHbond* firstHbond;

  double k;
  const double L;
  const double theta2;
  const double theta1;
  std::vector<char> isDonorH;
  std::vector<char> isAcceptor;
  std::vector<char> isBB;
  VerletCutoffsHAPairs* verletCutoffs;
  VerletCollectorHAPairs* verletCollector;

  typedef VL<VerletCutoffsHAPairs,VerletCollectorHAPairs> VerletList;
  VerletList* verlet;

  int nAdded_;
  int nRemoved_;
  bool doAngleConstraints;
  double Ecutoff;
  ForbidList* forbidlist;
  bool carefulChecking;

  void initializeConstraint( int h, int a, HBConstraint& hb );
  void initializeConstraint( int h, int a, BBHBConstraint& bbhb );
  void tightenConstraint( BBHBConstraint& bbhb );
  void tightenConstraint( HBConstraint& hb );

};

#endif /* HBMANAGER_H_ */
