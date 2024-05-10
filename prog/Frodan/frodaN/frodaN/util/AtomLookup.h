/*
 * AtomLookup.h
 *
 *  Created on: May 15, 2009
 *      Author: dwfarrel
 */

#ifndef ATOMLOOKUP_H_
#define ATOMLOOKUP_H_

class NeighborTable;
#include <vector>
#include <string>
#include <map>

class AtomLookup {
public:
  AtomLookup( const NeighborTable& nt_, const std::vector<std::string>& elem_ );
  ~AtomLookup() {}
  int operator()( const std::string& name ) { return lookup( name ); }
  void add( int index, const std::string& name );
  int lookup( const std::string& name );
  int neighH1( std::string name );
  void neighH3( std::string name, std::vector<int>& H );
  void clearfail() { failflag = false; }
  void clearatoms() { mapAtomNameToIndex.clear(); }
  bool fail() const { return failflag; }
private:
  const NeighborTable& nt;
  const std::vector<std::string>& elem;
  bool failflag;
  std::map<std::string,int> mapAtomNameToIndex;
};

#endif /* ATOMLOOKUP_H_ */
