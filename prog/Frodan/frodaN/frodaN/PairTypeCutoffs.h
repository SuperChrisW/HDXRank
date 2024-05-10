#ifndef PAIRTYPECUTOFFS_H_
#define PAIRTYPECUTOFFS_H_

#include <vector>

class PairTypeCutoffs
{
public:
  PairTypeCutoffs();
	virtual ~PairTypeCutoffs();

	void mapAtomsToTypes( const std::vector<int>& atomIndexToType_ );
  void setTypeRadii( const std::vector<double>& radii );
	void setPairTypeCutoff( int atomtype1, int atomtype2, double cutoff );

  double getCutoff( int p1, int p2 ) const {
    return cutoffLookup[atomIndexToType[p1]][atomIndexToType[p2]];
  }
	double getMaxCutoff() const { return maxCutoff; }

	int getAtomType( int atom ) const { return atomIndexToType[atom]; }

	void check();
private:
  int nAtoms;
  int nAtomTypes;
  std::vector<double> cutoffLookupStorage;
  std::vector<double*> cutoffLookup;
  std::vector<int> atomIndexToType;
  double maxCutoff;
};

#endif /*PAIRTYPECUTOFFS_H_*/
