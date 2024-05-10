/*
 * TargetMap.h
 *
 *  Created on: Aug 28, 2009
 *      Author: dwfarrel
 */

#ifndef TARGETMAP_H_
#define TARGETMAP_H_

#include <map>

class TargetMap {
public:
  TargetMap() {}
  virtual ~TargetMap() {}

  void insert( int isrc, int itarg ) {
    src2targ_[isrc] = itarg;
    targ2src_[itarg] = isrc;
  }
  void conv1To2( int isrc, int& itarg, bool& success ) const {
    std::map<int,int>::const_iterator itmap;
    itmap = src2targ_.find( isrc );
    success = ( itmap != src2targ_.end() );
    itarg = success ? itmap->second : -1;
  }
  void conv2To1( int itarg, int& isrc, bool& success ) const {
    std::map<int,int>::const_iterator itmap;
    itmap = targ2src_.find( itarg );
    success = ( itmap != targ2src_.end() );
    isrc = success ? itmap->second : -1;
  }

  const std::map<int,int>& src2targ() const { return src2targ_; }
  const std::map<int,int>& targ2src() const { return targ2src_; }

private:
  std::map<int,int> src2targ_;
  std::map<int,int> targ2src_;
};

#endif /* TARGETMAP_H_ */
