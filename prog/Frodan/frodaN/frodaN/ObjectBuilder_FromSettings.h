#ifndef OBJECTBUILDER_FROMSETTINGS_H_
#define OBJECTBUILDER_FROMSETTINGS_H_

class Settings;
class RigidUnitSystem;
class PDB;
class AmberPrmtop;
class ProteinInfo;
class NeighborTable;
class TargetMap;
#include <vector>
#include "Vec3.h"

class ObjectBuilder_FromSettings
{
public:
	ObjectBuilder_FromSettings( const Settings& settings );
	virtual ~ObjectBuilder_FromSettings();

  ProteinInfo* getProteinInfo();
  NeighborTable* getCovalentNeighborTable();
  RigidUnitSystem* getRigidUnitSystem();
  ProteinInfo* getTargetProteinInfo();
  NeighborTable* getTargetCovalentNeighborTable();
  TargetMap* getTargetMap();
  void getTargetCoords( std::vector<Vec3>& targCoords );

private:
  const Settings& settings;
  PDB* pdb;
  AmberPrmtop* prmtop;
  PDB* targpdb;
  AmberPrmtop* targprmtop;
  void initializeInputData();
  void initializeTargetData();

};

#endif /*OBJECTBUILDER_FROMSETTINGS_H_*/
