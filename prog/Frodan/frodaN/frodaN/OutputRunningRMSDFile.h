/*
 * OutputRunningRMSDFile.h
 *
 *  Created on: Jan 16, 2009
 *      Author: dwfarrel
 */

#ifndef OUTPUTRUNNINGRMSDFILE_H_
#define OUTPUTRUNNINGRMSDFILE_H_

class PDB;
class RigidUnitSystem;
class ProteinInfo;
#include <vector>
#include <string>
#include <fstream>
#include "Vec3.h"

class OutputRunningRMSDFile {
public:
  OutputRunningRMSDFile(
    const ProteinInfo* prot,
    const RigidUnitSystem* sys_,
    std::string fileprefix_ );
  virtual ~OutputRunningRMSDFile();
  void recordRMSD();
  void setOutputPeriod( int n ) { outputPeriod = n; }
private:
  const RigidUnitSystem* rigidUnitSystem;
  std::vector<int> caList;
  std::vector<int> resiList;
  std::vector<Vec3> initialPositions;
  std::ofstream rmsdfile;
  int outputPeriod;
  int iteration;
  double sum_of_msd;
  std::vector<double> caSumSquareDeviation;
  std::string fileprefix;
};

#endif /* OUTPUTRUNNINGRMSDFILE_H_ */
