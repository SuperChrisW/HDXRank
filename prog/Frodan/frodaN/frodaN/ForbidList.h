/*
 * ForbidList.h
 *
 *  Created on: Jun 21, 2010
 *      Author: dwfarrel
 */

#ifndef FORBIDLIST_H_
#define FORBIDLIST_H_

#include <map>
class Targeter;
class ProteinInfo;
class NeighborTable;

class ForbidList {
public:
  ForbidList( 
    int forbidTimeDuration_, 
    Targeter* targeter_,
    const ProteinInfo* prot_,
    const NeighborTable* nt_ ) :
    prot( prot_ ),
    nt( nt_ )
  {
    initialize( forbidTimeDuration_, targeter_ );
  }
  virtual ~ForbidList() {}
  void initialize( int forbidTimeDuration_, Targeter* targeter_ ) {
    forbidTimeDuration = forbidTimeDuration_;
    targeter = targeter_;
    timestampLookup.clear();
  }
  void insert( int p1, int p2 );
  void inserthb( int d, int a );

  bool checkforbid( int p1, int p2 ) const;
  void tidy();
private:
  const ProteinInfo* prot;
  const NeighborTable* nt;
  std::map< std::pair<int,int>, int > timestampLookup;
  int forbidTimeDuration;
  Targeter* targeter;
};

#endif /* FORBIDLIST_H_ */
