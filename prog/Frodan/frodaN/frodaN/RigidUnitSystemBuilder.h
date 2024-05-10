#ifndef RIGIDUNITSYSTEMBUILDER_H_
#define RIGIDUNITSYSTEMBUILDER_H_

class RigidUnitSystem;
#include <string>

class RigidUnitSystemBuilder
{
public:
	RigidUnitSystemBuilder() :
	  FIRSTBondFileInterpretation("pdb"),
	  hbondEnergyCutoff(0)
	{}
	virtual ~RigidUnitSystemBuilder() {}
	void setpdbfilename( std::string filename ) { pdbfilename = filename; }
	void setprmtopfilename( std::string filename ) { prmtopfilename = filename; }
	void setfirstrcfilename( std::string filename ) { firstrcfilename = filename; }
	void setfirstcovfilename( std::string filename ) { firstcovfilename = filename; }
	void setfirsthbondsfilename( std::string filename, double hbondEnergyCutoff_ ) { 
	  firsthbondsfilename = filename;
	  hbondEnergyCutoff = hbondEnergyCutoff_;
	}
	void setFIRSTBondFileInterpretation( std::string s ) {
	  FIRSTBondFileInterpretation = s;
	}
	void setrestartpdbfilename( std::string filename ) { restartpdbfilename = filename; }
	RigidUnitSystem *build();
private:
	std::string pdbfilename;
	std::string prmtopfilename;
	std::string firstrcfilename;
	std::string firstcovfilename;
	std::string restartpdbfilename;
	std::string FIRSTBondFileInterpretation;
	std::string firsthbondsfilename;
	double hbondEnergyCutoff;
};

#endif /*RIGIDUNITSYSTEMBUILDER_H_*/
