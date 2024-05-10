/*
 * PHManager.h
 *
 *  Created on: Aug 10, 2009
 *      Author: dwfarrel
 */



//Notes - for nucleotides:
//even within same nuc is ok.
//exclude pairs that have only 3 bonds between them or less.
//distance for nucs is 3.55
//need to identify bases at the beginning (isBase) vector char 0 or 1
//the base id is the residue id.
//need to be able to look up particular pair of bases and the corresponding pair of atoms



#ifndef PHMANAGER_H_
#define PHMANAGER_H_

#include "VL.h"
#include "NeighborTable.h"
#include "ProteinInfo.h"
#include <vector>

class HydrophobicPair {
public:
  HydrophobicPair() {}
  HydrophobicPair( int p1_, int p2_ ) : p1(p1_), p2(p2_) {}
  ~HydrophobicPair() {}
  int p1;
  int p2;
};

class VerletCollectorHydrophobicPairs {
public:
  VerletCollectorHydrophobicPairs() {}
  ~VerletCollectorHydrophobicPairs() {}
  void clearPairs() { pairs.clear(); }
  void addPair( int p1, int p2, double cutoff ) {
    pairs.push_back( HydrophobicPair(p1,p2) );
  }
  typedef std::vector< HydrophobicPair >::const_iterator const_iterator;
  const_iterator begin() const { return pairs.begin(); }
  const_iterator end() const { return pairs.end(); }
private:
  std::vector< HydrophobicPair > pairs;
};

class VerletCutoffsHydrophobicPairs {
public:
  VerletCutoffsHydrophobicPairs( const ProteinInfo& prot, const NeighborTable* nt_, const std::vector<char>* isNuc );
  ~VerletCutoffsHydrophobicPairs();
  void getCutoff( int p1, int p2, bool& isPairExcluded, double& cutoff ) const;
  double getMaxCutoff() const;
private:
  const NeighborTable* nt;
  const std::vector<char>* isNuc; //nucleic acids handling
  std::vector<int> lookupResIndexFromAtom;
};


//Added to handle nucleic acids following
//work of Holger Gohlke and Simone Fuller
class NucBasePair {
public:
  NucBasePair() {}

  NucBasePair( int b0, int b1 ) {
    base[0] = b0;
    base[1] = b1;
  }

  ~NucBasePair() {}

  bool operator<( const NucBasePair& other ) const {
    return std::lexicographical_compare( &base[0], &base[2], &other.base[0], &other.base[2] );
  }

  NucBasePair& operator=( const NucBasePair& other ) {
    base[0] = other.base[0];
    base[1] = other.base[1];
    return *this;
  }

  int base[2];
};

#include "Vec3.h"
#include "PHConstraint.h"
#include <map>
class ForbidList;

class PHManager
{
public:
  PHManager(
    const ProteinInfo& proteinInfo,
    const NeighborTable& nt,
    const std::vector<Vec3>& coords,
    PHContainer& phobes );
  virtual ~PHManager();

  void enableCarefulChecking() { carefulChecking = true; }
  void disableCarefulChecking() { carefulChecking = false; }
  void addConstraint( int p1, int p2, bool& success, bool breakable = true );
  void tighten();
  void findnew( bool breakable = true );
  void resetCounts() { nAdded_ = nRemoved_ = 0; }
  int nAdded() const { return nAdded_; }
  int nRemoved() const { return nRemoved_; }
  int nTotal() const { return phobes.size(); }
  bool isPairExcluded( int p1, int p2 ) const;

  void breakAllBreakable();
  void attachForbidList( ForbidList* forbidlist_ ) { forbidlist = forbidlist_; }

private:
  const ProteinInfo& prot;
  const NeighborTable& nt;
  const std::vector<Vec3>& coords;
  PHContainer& phobes;
  bool carefulChecking;

  //////
  //these variables added to handle nucleic acids following
  //work of Holger Gohlke and Simone Fuller
  std::vector<char> isNucAtom;
  std::vector<char> isNucBaseAtom;
  std::map<NucBasePair,HydrophobicPair> mapBasePairToPhobicPair;
  double Lnuc;
  //////

  std::vector<char> isAtomInHydrophobicVerletList;

  double k;
  const double L;
  VerletCutoffsHydrophobicPairs* verletCutoffs;
  VerletCollectorHydrophobicPairs* verletCollector;

  typedef VL<VerletCutoffsHydrophobicPairs,VerletCollectorHydrophobicPairs> VerletList;
  VerletList* verlet;

  int nAdded_;
  int nRemoved_;
  ForbidList* forbidlist;

  void initializeConstraint( int h, int a, PHConstraint& ph );
  void tightenConstraint( PHConstraint& hydrophobicConstraint );
  bool hasOnlyCHneighbors( int atom );

};

#endif /* PHMANAGER_H_ */
