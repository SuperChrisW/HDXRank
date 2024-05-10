#ifndef FRODACUTOFFSBUILDER_H_
#define FRODACUTOFFSBUILDER_H_

class PairTypeCutoffs;
class ProteinInfo;
class NeighborTable;
#include <vector>
#include <map>
#include <string>

class FrodaCutoffsBuilder
{
public:
  FrodaCutoffsBuilder(
	  const ProteinInfo &proteinInfo,
	  const NeighborTable &nt,
    double mismatchTol_,
    double scaleFactor );
	virtual ~FrodaCutoffsBuilder();
  PairTypeCutoffs* getPairTypeCutoffs() { return pairTypeCutoffs; }
private:
  PairTypeCutoffs* pairTypeCutoffs;
  bool polarGeometry;
  double polarHRadius;
  double vdwTol;
  double mismatchTol;
  std::vector<std::string> nameOfType;
  std::vector<double> radiusOfType;
  std::vector<signed char> chargeOfType;
  std::vector<bool> isPotentialDonorType;
  int nTypes;
  int nPairs;
  std::map<std::string, int> mapNameToType;
  std::vector<bool> isHydrogen;
  void setupTypes();
  void setupInteractionCutoffLookupTable();
  void assignTypes( const ProteinInfo &proteinInfo, const NeighborTable &nt );
  double calcInteractionCutoffForTypePair( int t1, int t2 ) const;
};

#endif /*FRODACUTOFFSBUILDER_H_*/
