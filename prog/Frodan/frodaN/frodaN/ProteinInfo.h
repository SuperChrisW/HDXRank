/*
 * ProteinInfo.h
 *
 *  Created on: Jul 22, 2009
 *      Author: dwfarrel
 */

#ifndef PROTEININFO_H_
#define PROTEININFO_H_

//class ProteinInfo contains the entire structure
// all atoms from 0..N-1.
// all residues from 0..N-1
//covalent bonding relationships are not implied, meaning that
//a residue may not be bonded to the next residue.  To find
//out, user must consult the neighbor table object

#include <vector>
#include <string>

class Atom;
class Resi;
class ProteinInfo;

class Atom {
public:
  Atom() {}
  Atom( int index__, std::string name__, std::string elem__, std::vector<Resi>* resi__, int resindex__ ) :
    index_(index__), name_(name__), elem_(elem__), resi_(resi__), resindex_(resindex__) {}
  ~Atom() {}
  void set( int index__, std::string name__, std::string elem__, std::vector<Resi>* resi__, int resindex__ ) {
    index_ = index__; name_ = name__; elem_ = elem__; resi_ = resi__; resindex_ = resindex__;
  }
  int index() const { return index_; }
  std::string name() const { return name_; }
  std::string elem() const { return elem_; }

  const Resi& resi() const { return (*resi_)[resindex_]; }
private:
  int index_;
  std::string name_;
  std::string elem_;
  std::vector<Resi>* resi_;
  int resindex_;
};

class Resi {
public:
  friend class ProteinInfo;
  Resi() {}
  ~Resi() {}
  Resi( int index__, std::string name__, std::vector<Atom>* atoms__, int beginatom__, int endatom__ ) :
    index_(index__), name_(name__), atoms_(atoms__), beginatom_(beginatom__), endatom_(endatom__) {}
  void set( int index__, std::string name__, std::vector<Atom>* atoms__, int beginatom__, int endatom__ ) {
    index_ = index__; name_ = name__; atoms_ = atoms__; beginatom_ = beginatom__; endatom_ = endatom__;
  }
  int natoms() const { return endatom_ - beginatom_; }
  int index() const { return index_; }
  std::string name() const { return name_; }
  typedef std::vector<Atom>::const_iterator const_iterator;
  const_iterator begin() const { return atoms_->begin() + beginatom_; }
  const_iterator end() const { return atoms_->begin() + endatom_; }
  const_iterator find( const std::string& atomname ) const {
    const_iterator it;
    for ( it = begin(); it != end(); it++ ) {
      if ( it->name() == atomname ) return it;
    }
    return it;
  }

private:
  int index_;
  std::string name_;
  std::vector<Atom>* atoms_;
  int beginatom_;
  int endatom_;
};

class PDB;
class AmberPrmtop;

class ProteinInfo {
public:
  ProteinInfo() {}
  ProteinInfo( const PDB& pdb );
  ProteinInfo( const AmberPrmtop& prmtop );
  virtual ~ProteinInfo() {}

  void startNewResidue( std::string name ) {
    int nextresindex = resi_.size();
    int beginatom = atom_.size();
    int endatom = beginatom;
    //if the last residue has no atoms, then this "new" residue overwrites
    //the old one.
    if ( nextresindex > 0 && resi_[nextresindex-1].natoms() == 0 ) {
      resi_[nextresindex-1] = Resi( nextresindex, name, &atom_, beginatom, endatom );
    }
    else {
      resi_.push_back( Resi( nextresindex, name, &atom_, beginatom, endatom ) );
    }
  }
  void addAtomToCurrentResidue( std::string name, std::string elem ) {
    int nextatomindex = atom_.size();
    int resindex = resi_.size() - 1;
    if ( resindex >= 0 ) {
      atom_.push_back( Atom( nextatomindex, name, elem, &resi_, resindex ) );
      resi_[resindex].endatom_++;
    }
  }

  const Resi& resi( int r ) const { return resi_[r]; }
  const Atom& atom( int a ) const { return atom_[a]; }

  typedef std::vector<Atom> Atoms;
  typedef std::vector<Resi> Residues;
  const Residues& residues() const { return resi_; }
  const Atoms& atoms() const { return atom_; }

  int natoms() const { return atom_.size(); }
  int nresi() const { return resi_.size(); }
private:
  std::vector<Atom> atom_;
  std::vector<Resi> resi_;

};

#endif /* PROTEININFO_H_ */
