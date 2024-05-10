/*
 * PDBTrajectoryInput.h
 *
 *  Created on: Jul 22, 2009
 *      Author: dwfarrel
 */

#ifndef PDBTRAJECTORYINPUT_H_
#define PDBTRAJECTORYINPUT_H_

#include "TrajectoryInput.h"
#include "Vec3.h"
#include "PDB.h"
#include <vector>
#include <string>

class PDBTrajectoryInput : public TrajectoryInput {
public:
  PDBTrajectoryInput();
  PDBTrajectoryInput( std::vector<std::string>& filenames );
  virtual ~PDBTrajectoryInput();

  void setFrameFilenames( std::vector<std::string>& filenames );
  void readNextFrame( std::vector<Vec3> &coords, bool &good );
  bool readNextFrame( std::vector<Vec3> &coords );
  void rewind();
private:
  bool readyToRead;
  bool verbose;
  int frame;
  std::vector<std::string> filenames;
  PDB pdb;
};

#endif /* PDBTRAJECTORYINPUT_H_ */
