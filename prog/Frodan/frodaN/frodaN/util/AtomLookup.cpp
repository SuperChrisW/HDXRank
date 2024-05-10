/*
 * AtomLookup.cpp
 *
 *  Created on: May 15, 2009
 *      Author: dwfarrel
 */

#include "AtomLookup.h"
#include "NeighborTable.h"
#include <cstdlib>
#include <iostream>

using namespace std;

AtomLookup::AtomLookup( const NeighborTable& nt_, const vector<string>& elem_ ) :
    nt( nt_ ),
    elem( elem_ )
{
      if ( nt.size() != elem.size() ) {
        cout << "Error, sizes not equal" << endl;
        exit(0);
      }
}

void AtomLookup::add( int index, const string& name ) {
    map<string,int>::iterator it;
    it = mapAtomNameToIndex.find( name );
    if ( it != mapAtomNameToIndex.end() ) {
      failflag = true;
      cout << "Error, cannot add duplicate atom name to lookup table " << index << " " << name << endl;
      exit(0);
    }
    else mapAtomNameToIndex[name] = index;
}

int AtomLookup::lookup( const string& name ) {
    map<string,int>::const_iterator it;
    it = mapAtomNameToIndex.find( name );
    if ( it == mapAtomNameToIndex.end() ) {
      failflag = true;
      return -1;
    }
    else return it->second;
}

int AtomLookup::neighH1( string name ) {
    int atom = lookup( name );
    if ( atom < 0 ) {
      failflag = true;
      return -1;
    }
    int H;
    int nNeigh = nt[atom].size();
    int count = 0;
    for ( int i = 0; i < nNeigh; i++ ) {
      int neigh = nt[atom][i];
      if ( elem[neigh] == "H" ) {
        count++;
        H = neigh;
      }
    }
    if ( count == 1 ) {
      return H;
    }
    else {
      failflag = true;
      return -1;
    }
}

void AtomLookup::neighH3( string name, vector<int>& H ) {
    H.resize(3);
    int atom = lookup( name );
    if ( atom < 0 ) {
      failflag = true;
      return;
    }
    int nNeigh = nt[atom].size();
    int count = 0;
    for ( int i = 0; i < nNeigh; i++ ) {
      int neigh = nt[atom][i];
      if ( elem[neigh] == "H" ) {
        if ( count < 3 ) H[count++] = neigh;
        else count++;
      }
    }
    if ( count != 3 ) failflag = true;
}
