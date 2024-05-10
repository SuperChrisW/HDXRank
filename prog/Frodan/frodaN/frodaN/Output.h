#ifndef OUTPUT_H_
#define OUTPUT_H_

class RigidUnitSystem;
class PDB;
class AmberTrajectory;
class RMSD;
#include "Observable.h"
#include <string>

class OutputSettings {
public:
  bool outputRMSDFiles;
  bool amberOutput;
  int newTrajPeriod;
  int outputConformerPeriod;
  double outputConfAtRMSDFromLast;
  std::string pdbfilename;
  std::string prmtopfilename;
};

class Output : public Observable
{
public:
  Output(
      const OutputSettings& settings,
      RigidUnitSystem *rigidUnitSystem );

  virtual ~Output();

  void notifyStructureReady();
private:
  const RigidUnitSystem* rigidUnitSystem;
  PDB* pdb;
  AmberTrajectory* traj;
  unsigned int iter;
  unsigned int snapshot;
  unsigned int outputFrequency;
  bool doRMSDFromLast;
  double rmsdTrigger;
  RMSD* rmsdFromLast;
  std::string pdbPrefix;

  void writeStructures();
  void setupPDBOutput( std::string inputpdbfilename );
  void setupAmberTrajOutput( std::string prmtopfilename );
  void writePDB();
  void writeAmberTraj();

};

#endif /*OUTPUT_H_*/
