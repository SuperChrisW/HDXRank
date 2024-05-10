#ifndef AMBERTRAJECTORYINPUT_H_
#define AMBERTRAJECTORYINPUT_H_

#include <fstream>
#include <string>
#include <vector>
#include "Vec3.h"
#include "AmberPrmtop.h"
#include "TrajectoryInput.h"

class AmberTrajectoryInput : public TrajectoryInput
{
public:
	AmberTrajectoryInput();
  AmberTrajectoryInput( const std::string &filename, const AmberPrmtop &prmtop );
  AmberTrajectoryInput( const std::string &filename, int natoms_, bool ifbox_ = false );
	virtual ~AmberTrajectoryInput();

	void open( const std::string &filename, const AmberPrmtop &prmtop );
	void open( const std::string &filename, int natoms_, bool ifbox_ = false );
  bool readNextFrame( std::vector<Vec3> &coords );
  void readNextFrame( std::vector<Vec3> &coords, bool &good );
  void readNextFrame( std::vector<Vec3> &coords, double &box0, double &box1, double &box2, bool &good );
  void close();
  void rewind();
private:
  std::ifstream trajfile;
  int nAtoms;
  int ifbox;
  int frame;
  bool verbose;
  std::string titlestring;

  void open( const std::string &filename );
};

#endif /*AMBERTRAJECTORYINPUT_H_*/
