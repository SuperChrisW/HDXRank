#include "Mismatch.h"
#include "Vect.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include "RigidUnitSystem.h"

using namespace std;

Mismatch::Mismatch() : 
  maxMismatch(0.0),
  isUpdated( false )
{
}

Mismatch::~Mismatch()
{
}

void Mismatch::receiveNotification( Observable *obs )
{
  const RigidUnitSystem *rigidUnitSystem = dynamic_cast<const RigidUnitSystem *>( obs );
  if ( rigidUnitSystem && rigidUnitSystem->AbsolutePositionsChanged() ) isUpdated = false;
  else if ( !rigidUnitSystem ) isUpdated = false;
}

void Mismatch::update() {  
  vector<double> mismatchList;  
  mismatchList.clear();
  for ( std::list<MismatchTerm*>::iterator it = mismatchTerms.begin();
        it != mismatchTerms.end();
        it++ ) {
    // 'it' is an iterator to a pointer
    // *it is a pointer
    mismatchList.push_back( (*it)->mismatch() );  
  }
  maxMismatch = *max_element( mismatchList.begin(), mismatchList.end() );
  isUpdated = true;
}

