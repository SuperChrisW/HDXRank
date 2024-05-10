#ifndef OUTPUTRIGIDUNITS_H_
#define OUTPUTRIGIDUNITS_H_

/*
class PDB;
#include <vector>
#include "Vec3.h"
#include <string>
#include "Observable.h"
class RigidUnitSystem;
class PerturbRelaxCycle;

class OutputRigidUnits : public Observer
{
public:
  OutputRigidUnits(
      std::string inputpdbfilename,
      const RigidUnitSystem *rigidUnitSystem_,
      const PerturbRelaxCycle *cycle_ );
  virtual ~OutputRigidUnits();

  void beforePerturb();
  void afterPerturb();
  void afterMinCycle();
  void finishedConformer();
  void write( std::string filename );
  void receiveNotification( Observable *obs );
private:
  PDB *pdb;
  const RigidUnitSystem *rigidUnitSystem;
  const PerturbRelaxCycle *cycle;
  unsigned int countConf;
  unsigned int minCycleCount;
};
*/
#endif /*OUTPUTRIGIDUNITS_H_*/
