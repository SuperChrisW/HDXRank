/*
 * CellList.cpp
 *
 *  Created on: Feb 10, 2010
 *      Author: dwfarrel
 */

#include "CellList.h"
#include <cmath>

using namespace std;

CellList::CellList( double dist_ ) {
  setMaxDist( dist_ );
  neighCells.reserve( 5*5*5 );
  initializePairFinder();
}

CellList::~CellList() {
}

void CellList::update() {
  //empty the cells, but keep the allocated space of the cells
  for ( map< CellID, Cell >::iterator it = cells.begin();
        it != cells.end(); it++ ) {
    it->second.clear();
  }

  //fill cells
  size_t N = points.size();
  for ( size_t i = 0; i < N; i++ ) {
    const Vec3* r = points[i].position;
    CellID cellid;
    cellid.bin[0] = bin(r->x);
    cellid.bin[1] = bin(r->y);
    cellid.bin[2] = bin(r->z);
    cells[cellid].push_back( points[i].id );
  }

  //prune empty cells
  //This is vital, because the other routines rely on the fact that
  //empty cells are removed from the map.
  map< CellID, Cell >::iterator it = cells.begin();
  map< CellID, Cell >::iterator temp;
  while ( it != cells.end() ) {
    if ( it->second.empty() ) {
      temp = it;
      temp++;
      cells.erase( it );
      it = temp;
    }
    else it++;
  }

  initializePairFinder();
}

void CellList::setMaxDist( double dist_ ) {
  edgelength = dist_/2.0;
}

void CellList::insert( int id, const Vec3* position ) {
  points.push_back( Point( id, position ) );
}

void CellList::initializePairFinder() {
  firsttime = true;
  i_neighCellPoint = 0;
  i_neighCell = 0;
  i_refCellPoint = 0;
  ref_it = cells.begin();

  //It is possible that the list is totally empty,
  //the ref_it is already at the end of the cell list.
  if ( ref_it != cells.end() ) {
    refCell = &ref_it->second;
    setNeighborCells();
  }
  else {
    refCell = NULL;
    neighCells.clear();
  }

}

void CellList::setNeighborCells() {

  neighCells.clear();

  CellID cellid = ref_it->first;
  int binxstart = cellid.bin[0] - 2;
  int binxstop = cellid.bin[0] + 2;
  int binystart = cellid.bin[1] - 2;
  int binystop = cellid.bin[1] + 2;
  int binzstart = cellid.bin[2] - 2;
  int binzstop = cellid.bin[2] + 2;


  for ( int binx = binxstart; binx <= binxstop; binx++ ) {
    for ( int biny = binystart; biny <= binystop; biny++ ) {
      for ( int binz = binzstart; binz <= binzstop; binz++ ) {
        map< CellID, Cell >::const_iterator it = cells.find( CellID( binx, biny, binz) );
        if ( it == cells.end() ) continue;
        neighCells.push_back( &it->second );
      }
    }
  }
}

void CellList::getNextPair( int& id1, int& id2, bool& valid ) {
  //valid is initially set to false
  valid = false;

  //if we're at the end, then there are no more pairs. Trying to
  //advance further will cause errors, So just return.
  if ( ref_it == cells.end() ) return;

  //Advance to next valid pair.
  //But if it's the first time that a pair is being requested,
  //don't advance.
  if ( !firsttime ) {
    if ( ++i_neighCellPoint >= (*neighCells[i_neighCell]).size() ) {
      i_neighCellPoint = 0;
      if ( ++i_neighCell >= neighCells.size() ) {
        i_neighCell = 0;
        if ( ++i_refCellPoint >= (*refCell).size() ) {
          i_refCellPoint = 0;
          if ( ++ref_it != cells.end() ) {
            refCell = &ref_it->second;
            setNeighborCells();
          }
        }
      }
    }
  }
  else {
    firsttime = false;
  }

  //Only output a valid pair if the advance did not take us
  //to the end.
  if ( ref_it != cells.end() ) {
    id1 = (*refCell)[i_refCellPoint];
    id2 = (*neighCells[i_neighCell])[i_neighCellPoint];
    valid = true;
  }

}

//probably need to verify that the bin calculation does not cause an integer overflow,
//for example if the coord/length is a number that is larger than maxint or smaller than
//negative maxint.
//But for now we assume that this will not happen.
int CellList::bin( const double& coord ) const {
  return static_cast<int>( floor( coord/edgelength ) );
}
