#ifndef PERTURBRELAXCYCLE_H_
#define PERTURBRELAXCYCLE_H_

#include "PerturbRelaxCycleBase.h"
class MinimizeSystem;
class GenericMap;
class RandomRotorPerturber;
class RandomCenterPerturber;
class SymmetricPerturber;
class GlobalMotionRemover;
class MomentumPerturber;
class PhaseSpacePathLengthIntegrator;
class RMSD;
class RigidUnitSystem;
class ConstraintEnforcingPotential;
class Settings;
class EMStructure;
class CorrPerturber;
class MapPerturber;
class NeighborTable;
class ProteinInfo;
class Output;

class PerturbRelaxCycle : public PerturbRelaxCycleBase
{
public:
	PerturbRelaxCycle( const Settings& settings );
	virtual ~PerturbRelaxCycle();

	RigidUnitSystem* sys;
	ProteinInfo* prot;
	NeighborTable* nt;
	ConstraintEnforcingPotential* cep;
	MinimizeSystem *minim;
	GenericMap *map;
	RandomRotorPerturber *randomRotorPerturber;
	RandomCenterPerturber *randomCenterPerturber;
	SymmetricPerturber *symPert;
	GlobalMotionRemover *globalMotionRemover;
	MomentumPerturber *mom;
	PhaseSpacePathLengthIntegrator *pathLengthIntegrator;
	RMSD *rmsdFromInitial;
	EMStructure *myEM;
	CorrPerturber *corrPert;
	MapPerturber *mapPert;
	Output* output;

	void outputSummaryLine();
};

#endif /*PERTURBRELAXCYCLE_H_*/
