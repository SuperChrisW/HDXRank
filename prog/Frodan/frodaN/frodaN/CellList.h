/*
 * CellList.h
 *
 *  Created on: Feb 10, 2010
 *      Author: dwfarrel
 */

#ifndef CELLLIST_H_
#define CELLLIST_H_

#include <map>
#include <vector>
#include "Vec3.h"

class CellList {
public:
  class Point;
  class CellID;
  typedef std::vector<int> Cell;

  CellList( double dist );
  virtual ~CellList();

  void setMaxDist( double dist );
  void insert( int id, const Vec3* position );
  void update();
  void getNextPair( int& id1, int& id2, bool& valid );
private:
  std::vector<Point> points;
  std::map< CellID, Cell > cells;
  double edgelength;
  bool nextPairReady;
  std::map< CellID, Cell >::const_iterator ref_it;
  const Cell* refCell;
  std::vector<const Cell*> neighCells;
  size_t i_neighCell;
  size_t i_refCellPoint;
  size_t i_neighCellPoint;
  bool firsttime;

  int bin( const double& coord ) const;
  void setNeighborCells();
  void initializePairFinder();

};

class CellList::Point {
public:
  Point() {}
  Point( int id_, const Vec3* position_ ) :
    id( id_ ),
    position( position_ ) {}
  int id;
  const Vec3 *position;
};

#include <algorithm>

class CellList::CellID {
public:
  CellID() {}

  CellID( const int& x, const int& y, const int& z ) {
    bin[0] = x; bin[1] = y; bin[2] = z;
  }

  ~CellID() {}

  bool operator<( const CellID& other ) const {
    return std::lexicographical_compare( &bin[0], &bin[3], &other.bin[0], &other.bin[3] );
  }

  CellID& operator=( const CellID& other ) {
    bin[0] = other.bin[0];
    bin[1] = other.bin[1];
    bin[2] = other.bin[2];
    return *this;
  }

  int bin[3];
};

#endif /* CELLLIST_H_ */
