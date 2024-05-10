#ifndef AMBERTRAJECTORY_H_
#define AMBERTRAJECTORY_H_

#include "Vec3.h"
#include <string>
#include <vector>
class AmberPrmtop;

class AmberTrajectory
{
public:
  AmberTrajectory();
  virtual ~AmberTrajectory();

  void writeRestartFile( const std::string &filename_, bool box, const std::vector<Vec3> &coords );

  void initializeOutputFile( const std::string &filename_, bool box );
  void append( const std::vector<Vec3> &coords );

private:
  std::string filename;
  bool readyToAppend;
  bool ifbox;
  int width;
  int nCols;
  int precision;
};

#endif /*AMBERTRAJECTORY_H_*/
