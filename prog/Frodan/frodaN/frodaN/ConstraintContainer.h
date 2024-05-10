/*
 * ConstraintContainer.h
 *
 *  Created on: Aug 7, 2009
 *      Author: dwfarrel
 */

#ifndef CONSTRAINTCONTAINER_H_
#define CONSTRAINTCONTAINER_H_

#include <vector>
#include <map>
#include "Energy.h"
#include "Gradient.h"
#include "Observable.h"

class InsertDeleteEventData {
public:
  short int code;
  int p1;
  int p2;
};

class ConstraintContainerAbstract : public Observable {
public:
  virtual ~ConstraintContainerAbstract() {}
  virtual const InsertDeleteEventData& getInsertDeleteEventData() const = 0;
};

template <class T>
class ConstraintContainer :
  public EnergyTerm,
  public GradientTerm_P,
  public ConstraintContainerAbstract
{
public:
  ConstraintContainer( int natoms ) : lookupConstraintIndex( natoms ), enabled( true ) {}
  virtual ~ConstraintContainer() {}

  double energy();
  void addToGradient_P(
    std::vector<Vec3> &dV_dr_P,
    std::vector<SecondDerivative> &secondDerivative_P );

  typedef typename std::vector<T>::iterator iterator;
  typedef typename std::vector<T>::const_iterator const_iterator;

  iterator begin() { return constraints.begin(); }
  iterator end() { return constraints.end(); }
  const_iterator begin() const { return constraints.begin(); }
  const_iterator end() const { return constraints.end(); }
  size_t size() const { return constraints.size(); }

  iterator find( int p1, int p2 );
  const_iterator find( int p1, int p2 ) const;
  void lookupConstraintPartners( int p, std::vector<int>& partners ) const;
  void erase( int p1, int p2 );
  iterator erase( iterator it );
  iterator insert( int p1, int p2, const T& constraint );
  void clear() {
    constraints.clear();
    size_t natoms = lookupConstraintIndex.size();
    lookupConstraintIndex.clear();
    lookupConstraintIndex.resize( natoms );
    indexPairs.clear();
  }
  void enable() { enabled = true; }
  void disable() { enabled = false; }
  const InsertDeleteEventData& getInsertDeleteEventData() const { return notify; }
protected:
  std::vector< std::map< int, size_t > > lookupConstraintIndex;
  std::vector<T> constraints;
  std::vector< std::pair<int,int> > indexPairs;
  bool enabled;
  InsertDeleteEventData notify;
};

template <class T>
inline typename ConstraintContainer<T>::const_iterator ConstraintContainer<T>::find( int p1, int p2 ) const {
  if ( static_cast<size_t>( p1 ) >= lookupConstraintIndex.size() ||
       static_cast<size_t>( p2 ) >= lookupConstraintIndex.size() ) return constraints.end();
  const std::map< int, size_t >& lookupPartners = lookupConstraintIndex[p1];
  std::map< int, size_t >::const_iterator it = lookupPartners.find( p2 );
  if ( it == lookupPartners.end() ) { return constraints.end(); }
  return constraints.begin() + it->second;
}

template <class T>
inline typename ConstraintContainer<T>::iterator ConstraintContainer<T>::find( int p1, int p2 ) {
  if ( static_cast<size_t>( p1 ) >= lookupConstraintIndex.size() ||
       static_cast<size_t>( p2 ) >= lookupConstraintIndex.size() ) return constraints.end();
  std::map< int, size_t >& lookupPartners = lookupConstraintIndex[p1];
  std::map< int, size_t >::iterator it = lookupPartners.find( p2 );
  if ( it == lookupPartners.end() ) { return constraints.end(); }
  return constraints.begin() + it->second;
}

template <class T>
inline void ConstraintContainer<T>::lookupConstraintPartners(
  int p, std::vector<int>& partners ) const {
  partners.clear();
  if ( static_cast<size_t>( p ) >= lookupConstraintIndex.size() ) return;
  const std::map< int, size_t >& lookupPartners = lookupConstraintIndex[p];
  std::map< int, size_t >::const_iterator it;
  for ( it = lookupPartners.begin(); it != lookupPartners.end(); it++ ) {
    partners.push_back( it->first );
  }
}

template <class T>
inline void ConstraintContainer<T>::erase( int p1, int p2 ) {
  iterator it = find( p1, p2 );
  if ( it != end() ) erase( it );
}

template <class T>
inline typename ConstraintContainer<T>::iterator ConstraintContainer<T>::erase( iterator it ) {
  // this code will fail if p1 == p2.  For this reason, in the
  // "insert" member function, constraints with
  // p1 == p2 are rejected.

  size_t index = it - constraints.begin();

  //lookup the two atom indices for this constraint
  int p1 = indexPairs[index].first;
  int p2 = indexPairs[index].second;

  notify.code = -1;
  notify.p1 = p1;
  notify.p2 = p2;
  notifyObservers();

  //erase the entries in the lookup table for this constraint
  std::map< int, size_t >& partners1 = lookupConstraintIndex[p1];
  std::map< int, size_t >& partners2 = lookupConstraintIndex[p2];
  std::map< int, size_t >::iterator it1 = partners1.find( p2 );
  std::map< int, size_t >::iterator it2 = partners2.find( p1 );
  partners1.erase( it1 );
  partners2.erase( it2 );

  //move the last constraint into the place of the deleted constraint,
  //and shrink the constraint list by one.
  *it = constraints.back();
  indexPairs[index] = indexPairs.back();
  constraints.resize( constraints.size() - 1 );
  indexPairs.resize( indexPairs.size() - 1 );

  if ( index < constraints.size() ) {
    //since we have changed the index of a constraint (the last constraint has
    //moved to replace the erased constraint),
    //we need to locate its entries
    //in the lookup table and give them the new constraint index.
    p1 = indexPairs[index].first;
    p2 = indexPairs[index].second;
    lookupConstraintIndex[p1][p2] = lookupConstraintIndex[p2][p1] = index;
    return constraints.begin() + index;
  }
  else return constraints.end();
}

template <class T>
inline typename ConstraintContainer<T>::iterator ConstraintContainer<T>::insert( int p1, int p2, const T& constraint ) {
  // this code will fail if p1 == p2.  For this reason, we reject pairs with p1 == p2

  //if p1 or p2 is out of bounds, cannot insert.
  if ( static_cast<size_t>( p1 ) >= lookupConstraintIndex.size() ||
       static_cast<size_t>( p2 ) >= lookupConstraintIndex.size() ) return constraints.end();
  //also do not allow inserting of p1 == p2 (cannot have a constraint with itself)
  if ( p1 == p2 ) return constraints.end();

  notify.code = 1;
  notify.p1 = p1;
  notify.p2 = p2;
  notifyObservers();


  std::map< int, size_t >& partners1 = lookupConstraintIndex[p1];
  std::map< int, size_t >& partners2 = lookupConstraintIndex[p2];
  std::map< int, size_t >::iterator it1 = partners1.find( p2 );
  if ( it1 == partners1.end() ) {
    partners1.insert( std::pair<int, size_t>( p2, constraints.size() ) );
    partners2.insert( std::pair<int, size_t>( p1, constraints.size() ) );
    constraints.push_back( constraint );
    indexPairs.push_back( std::pair<int,int>( p1, p2) );
    return constraints.end() - 1;
  }
  else {
    constraints[it1->second] = constraint;
    return constraints.begin() + it1->second;
  }
}

template <class T>
inline double ConstraintContainer<T>::energy() {
  if ( !enabled ) return 0;
  double E = 0;
  typename std::vector<T>::iterator it;
  for ( it = constraints.begin(); it != constraints.end(); it++ ) {
    E += it->energy();
  }
  return E;
}

template <class T>
inline void ConstraintContainer<T>::addToGradient_P(
  std::vector<Vec3> &dV_dr_P,
  std::vector<SecondDerivative> &secondDerivative_P ) {
  if ( !enabled ) return;
  typename std::vector<T>::iterator it;
  for ( it = constraints.begin(); it != constraints.end(); it++ ) {
    it->addToGradient_P( dV_dr_P, secondDerivative_P );
  }
}

#endif /* CONSTRAINTCONTAINER_H_ */
