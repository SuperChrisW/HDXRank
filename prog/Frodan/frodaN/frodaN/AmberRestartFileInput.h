#ifndef AMBERRESTARTFILEINPUT_H_
#define AMBERRESTARTFILEINPUT_H_

#include <fstream>
#include <string>
#include <vector>
#include "Vec3.h"
#include "AmberPrmtop.h"

class AmberRestartFileInput
{
public:
  AmberRestartFileInput();
	virtual ~AmberRestartFileInput();

  void read( const std::string &filename, const AmberPrmtop &prmtop, std::vector<Vec3> &coords, bool &good );
};

#endif /*AMBERTRAJECTORYINPUT_H_*/
