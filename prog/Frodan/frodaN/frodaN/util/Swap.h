/*
 * Swap.h
 *
 *  Created on: May 15, 2009
 *      Author: dwfarrel
 */

#ifndef SWAP_H_
#define SWAP_H_

#include <vector>
#include <map>
#include "Vec3.h"
#include "Fit.h"

class SwapBase {
public:
  SwapBase() {}
  virtual ~SwapBase()=0;
  const std::vector<int>& getIndices() const { return indices; }
  const std::vector< std::vector<int> >& getAltMaps() const { return altmaps; }
protected:
  std::vector<int> indices;
  std::vector< std::vector<int> > altmaps;
};

class SwapSet_2Fold : public SwapBase {
public:
  SwapSet_2Fold();
  virtual ~SwapSet_2Fold() {}
  void clear();
  void addPair( int i, int j );
};

class SwapPoint_3Fold : public SwapBase {
public:
  SwapPoint_3Fold();
  virtual ~SwapPoint_3Fold() {}
  void setTriple( int i, int j, int k );
  void setTriple( const std::vector<int>& triple );
};

class CoordinateSwapper {
public:
  CoordinateSwapper(
    std::vector<Vec3>& state1coords_,
    const std::vector<Vec3>& state2coords_,
    const std::map<int,int>* atommap_ = NULL );

  ~CoordinateSwapper() {}

  void swap(
    const std::vector<int>& fittingSubset1,
    const SwapBase& swapBase );
private:
  std::vector<Vec3>& state1coords;
  const std::vector<Vec3>& state2coords;
  const std::map< int, int >* atommap;
  Fit fit;

  void conv1to2(
    const std::vector<int>& indices1,
    std::vector<int>& indices2,
    bool& success ) const ;
  void swapIndices(
    const std::vector<int>& indices,
    const std::vector<int>& remapping,
    std::vector<int>& swappedIndices ) const;

  void extractCoords(
    const std::vector<Vec3>& coords,
    const std::vector<int>& indices,
    std::vector<Vec3>& extractedCoords ) const;

  void setCoords(
    std::vector<Vec3>& coords,
    const std::vector<int>& indices,
    const std::vector<Vec3>& newCoords ) const;
  double energy( const std::vector<Vec3>& positions1, const std::vector<Vec3>& positions2 ) const;

  void pickBest( const std::vector<int>& indices,
    const std::vector< std::vector<int> >& altmaps,
    const std::vector<Vec3>& reference );

};

#endif /* SWAP_H_ */
