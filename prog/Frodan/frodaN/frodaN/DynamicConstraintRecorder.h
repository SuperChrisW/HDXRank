/*
 * DynamicConstraintRecorder.h
 *
 *  Created on: Jun 18, 2010
 *      Author: dwfarrel
 */

#ifndef DYNAMICCONSTRAINTRECORDER_H_
#define DYNAMICCONSTRAINTRECORDER_H_

#include "Observable.h"
#include <fstream>
#include "ConstraintContainer.h"

template < class T >
class DynamicConstraintRecorder : public Observer {
public:
  DynamicConstraintRecorder( std::string filename, ConstraintContainerAbstract* c_, const T* t_ ) :
    c( c_ ),
    t( t_ ) {

    o.open( filename.c_str(), std::ios::out );
    c->registerObserver( this );
  }
  virtual ~DynamicConstraintRecorder() {
    o.close();
  }
  void receiveNotification( Observable *observable ) {
    if ( observable == c ) {
      o << c->getInsertDeleteEventData().p1 << " " <<
            c->getInsertDeleteEventData().p2 << " " <<
            c->getInsertDeleteEventData().code << " " <<
            t->getIteration() << std::endl;
    }
  }
private:
  std::ofstream o;
  ConstraintContainerAbstract* c;
  const T* t;
};

#endif /* DYNAMICCONSTRAINTRECORDER_H_ */
