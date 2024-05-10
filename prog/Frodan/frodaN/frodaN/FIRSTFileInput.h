#ifndef FIRSTFILEINPUT_H_
#define FIRSTFILEINPUT_H_

class UniqueIDtoArrayIndexConverter;
class PDB;
class RigidUnitSystem;
class NeighborTable;
//class Ropes;
#include <string>
#include <vector>
#include "Vec3.h"

class FIRSTFileInput
{
public:
	FIRSTFileInput( std::string pdbfilename, std::string FIRSTBondFileInterpretation );
  FIRSTFileInput( const PDB* pdb, std::string FIRSTBondFileInterpretation );

	virtual ~FIRSTFileInput();

	//Ropes *buildRopes_FromFIRSThphobes(
	//    std::string phobesfilename,
	//    const RigidUnitSystem *rigidUnitSystem );
	NeighborTable* buildCovalentNeighborTable_FromFIRSTcov(
	    std::string firstcovfilename );
	NeighborTable *buildHbondNeighborTable( std::string filename, double energyCutoff );
  const std::vector<Vec3>& getPDBAtomPositions() { return initialPoints; }
  std::vector< std::vector<int> > *buildRUtoPlist_FromFIRSTdatafile(
      std::string filename, const NeighborTable& neighborTable );

private:
  PDB* internalpdb;
  const PDB* pdb;
  UniqueIDtoArrayIndexConverter *converter;
  std::vector<Vec3> initialPoints;
  void setup( std::string FIRSTBondFileInterpretation );

};

#endif /*FIRSTFILEINPUT_H_*/
