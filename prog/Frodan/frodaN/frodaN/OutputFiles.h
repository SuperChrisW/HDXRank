#ifndef OUTPUTFILES_H_
#define OUTPUTFILES_H_

class RigidUnitSystem;
class Targeter;
class ProteinInfo;
class OutputRunningRMSDFile;
#include "Observable.h"
#include <string>
#include <fstream>

class OutputFiles : public Observer
{
public:
  OutputFiles( const RigidUnitSystem *rigidUnitSystem_ );
  virtual ~OutputFiles();

  void setupRMSDToTargetOutput( const Targeter *targ_ );
  void setupRunningRMSD( const ProteinInfo *prot, std::string pdbfilename );

  void receiveNotification( Observable* obs );
  void writeFiles();
private:
  const RigidUnitSystem *rigidUnitSystem;
  unsigned int snapshot;
  const Targeter* targ;
  OutputRunningRMSDFile* outputRunningRMSDFile;
  std::string filenameRMSDToTarget;
  std::ofstream fileRMSDToTarget;

  void writeRMSDToTarget();

};

#endif /*OUTPUTFILES_H_*/
