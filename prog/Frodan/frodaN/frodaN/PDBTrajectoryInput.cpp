/*
 * PDBTrajectoryInput.cpp
 *
 *  Created on: Jul 22, 2009
 *      Author: dwfarrel
 */

#include "PDBTrajectoryInput.h"
#include <iostream>
#include <iomanip>

using namespace std;

PDBTrajectoryInput::PDBTrajectoryInput() :
  verbose( true ),
  frame(0)
{
}

PDBTrajectoryInput::~PDBTrajectoryInput() {
}

PDBTrajectoryInput::PDBTrajectoryInput( std::vector<std::string>& filenames_ ) :
  verbose( true ) {
  setFrameFilenames( filenames_ );
}

void PDBTrajectoryInput::setFrameFilenames( std::vector<std::string>& filenames_ ) {
  filenames = filenames_;
  frame = 0;
}

void PDBTrajectoryInput::rewind() {
  frame = 0;
}

bool PDBTrajectoryInput::readNextFrame( std::vector<Vec3> &coords ) {
  bool good;
  readNextFrame( coords, good );
  return good;
}

void PDBTrajectoryInput::readNextFrame( std::vector<Vec3> &coords, bool &good ) {

  int nFrames = filenames.size();
  if ( frame < nFrames ) {
    pdb.read( filenames[frame] );
    pdb.getPositions( coords );
    good = true;
  }
  else {
    good = false;
  }

  if ( verbose ) {
    if ( good ) {
      if ( (frame%50)==0 ) {
        cout << "\n" << setw(6) << right << frame << flush;
      }
      cout << "." << flush;
    }
    else {
      cout << "\n" << frame << " frames successfully read from trajectory." << endl;
    }
  }
  frame++;

}
