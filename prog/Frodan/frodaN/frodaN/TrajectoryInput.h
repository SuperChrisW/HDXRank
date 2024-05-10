/*
 * TrajectoryInput.h
 *
 *  Created on: Jul 22, 2009
 *      Author: dwfarrel
 */

#ifndef TRAJECTORYINPUT_H_
#define TRAJECTORYINPUT_H_

#include <vector>
#include "Vec3.h"

class TrajectoryInput {
public:
  TrajectoryInput() {}
  virtual ~TrajectoryInput() {}
  virtual bool readNextFrame( std::vector<Vec3> &coords ) = 0;
  virtual void readNextFrame( std::vector<Vec3> &coords, bool &good ) = 0;
  virtual void close() {}
  virtual void rewind() = 0;
};

#endif /* TRAJECTORYINPUT_H_ */
