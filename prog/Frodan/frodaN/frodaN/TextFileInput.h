#ifndef TEXTFILEINPUT_H_
#define TEXTFILEINPUT_H_

#include <string>
class Ropes;
class RigidUnitSystem;
class NeighborTable;
class PairTypeCutoffs;
class AmberPrmtop;
class TargetEnergy;
class TargetMap;
class HBManager;
class PHManager;
#include <map>
#include <vector>
#include <set>

class TextFileInput
{
public:
	TextFileInput();
	virtual ~TextFileInput();

	//Ropes *buildRopes(
	//    std::string filename,
	//    const RigidUnitSystem *rigidUnitSystem );

	NeighborTable *buildCovalentNeighborTable_FromPrmtop(
	    const AmberPrmtop& prmtop );

	void buildAtomTypes( std::string filename, std::vector<int>& atomtypes );

	PairTypeCutoffs *buildPairTypeCutoffs(
	    std::string atomTypesAssignmentFilename,
	    std::string pairTypeCutoffsFilename,
	    int nAtoms,
	    double scaleFactor = 1.0 );

	void readConstraintsFile(
	  std::string filename,
	  std::vector<int>& firstatom,
	  std::vector<int>& secondatom,
	  std::vector<double>& dist );
	void readPairsFile(
	  std::string filename,
	  std::vector<int>& firstatom,
	  std::vector<int>& secondatom,
          std::string format = "index0" );
        void readHBondFile_index0( std::string filename, HBManager *hbManager );
        void readHydrophobicsFile_index0( std::string filename, PHManager *phManager );


	TargetEnergy *buildTargetEnergy(
	    std::string filename,
	    RigidUnitSystem *rigidUnitSystem );

	std::map<int,int>* buildAtomMap( std::string mapFilename );
	TargetMap* buildTargetMap( std::string mapFilename );

	void readAtomSet( std::string filename, std::set<int>& atomset );

};

#endif /*TEXTFILEINPUT_H_*/
