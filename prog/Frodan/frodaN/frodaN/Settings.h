#ifndef SETTINGS_H_
#define SETTINGS_H_

#include <string>
#include "Output.h"
#define TIXML_USE_TICPP
#include "ticpp.h"

class Settings
{
public:
  Settings() {}
	Settings( int argc, char **argv );
	virtual ~Settings();

	std::string runtype;
  bool unbreakableInitialCon;
  int switchToBreak;
  int switchToAdd;
  int switchOff;
  int switchBackToIterZero;

  bool outputconstraintlists;

  bool targSubsetHeavy;
  double targDelta;
  std::string targ;
  std::string commonConstraintHandling;
  std::string noncommonConstraintHandling;

  double pairCutoffScaleFactor;
  bool doBacktrack;
  bool doAutoOverrideMinDist;
  bool doRedetermineConstraints;
  bool targDynamicConstraints;
  bool mergeConstraintsForTargeting;


  class Files {
   public:
    std::string prmtop;
    std::string firstcov;
    std::string pdb;
    std::string firsthbonds;
    std::string firstphobes;
    std::string firstrc;
    std::string restartpdb;
    double hbondEnergyCutoff;
    std::string FIRSTBondFileInterpretation;
    bool rigidHbonds;
    std::string symMatrices;
    std::string ezdMap;
    std::string lessThanConstraintsFile;
    std::string angleLessThanConstraintsFile;
    std::string overrideMinDistFile;
    std::string targOverrideMinDistFile;
    std::string atomTypesAssignment;
    std::string pairTypeCutoffs;
    std::string targfirstcov;
    std::string targatomtypes;
    std::string targfirstrc;
    std::string atommap;
    std::string targpdb;
    std::string targprmtop;
    std::string targAmberRestart;
    std::string targsubset;
    std::string breakableconstraints_atomlist;
    std::string hbonds_index0;
    std::string hydrophobics_index0;
  };

  class Energy {
  public:
    bool doOverlapEnergy;
    bool doSymmetryEnergy;
    bool doRama;
    bool doTorsion;
    bool doHbond;
    bool doHydrophobic;
    bool doHbondAngles;
  };

  class PerturbRelax {
   public:
    std::string tolType;
    double tol;
    int Nsteps;
    int NminimizationSteps;
    int startStep;
    double randomCenterPerturbationSize;
    bool doRandomCenterPerturbation;
    double randomRotorPerturbationSize;
    bool doRandomRotorPerturbation;
    bool doMomentumPerturbation;
    double momentumScaleFactor;
    bool removeGlobalMotion;
    bool stopAtPhysicalTimeLimit;
    double physicalTimeLimitPS;
    double symmetricPerturbationSize;
    bool doSymmetricPerturbation;
    double mapPerturbationSize;
    bool doMapPerturbation;
    double corrPerturbationSize;
    bool doCorrPerturbation;
    double emResolution;
    double emCutoff;
    std::string fitSelection;
  };


  Files files;
  Energy energy;
  PerturbRelax perturbRelax;
  OutputSettings output;
  unsigned long int seed;
  std::string repulsionType;
  std::string seedstring;

private:
  void sanityCheck() const;  // check options for consistency
  void parseXMLSettings_runfixedconstraints( ticpp::Element* e );
  void parseXMLSettings_runtarget( ticpp::Element* e );
  void parseXMLSettings_random( ticpp::Element* e );
  void parseXMLSettings_constraints( ticpp::Element* e );
  void parseXMLSettings_output( ticpp::Element* e );
  void parseXMLSettings_modelfiles( ticpp::Element* e );
  void parseXMLSettings_targetfiles( ticpp::Element* e );
  void parseXml( std::string filename );

};

#endif /*SETTINGS_H_*/
