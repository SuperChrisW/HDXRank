#ifndef ATOMCATEGORIES_H_
#define ATOMCATEGORIES_H_

#include <vector>
#include <string>
#include <map>
#include <set>
class ProteinInfo;
class NeighborTable;

class AtomCategories
{
public:
	AtomCategories( const NeighborTable *neighborTable_, const ProteinInfo *prot );
	virtual ~AtomCategories();
	int getAtomTypeForAtom( int atom ) const { return atomIndexToType[atom]; }

private:
  const NeighborTable *neighborTable;
  const ProteinInfo *prot;
  std::vector<std::string> types;
  std::map<std::string, int> mapTypeNameToTypeIndex;
  std::vector<int> atomIndexToType;

  std::set<std::string> proteinResidueNames;
  std::set<std::string> aromaticCarbonNames;
  std::set<std::string> amideNitrogenNames;
  std::set<std::string> carboxylOxygenNames;

  void defineTypes();
  void assignTypesToAtoms();
  void getProteinAtomType( int atom, int& atomtype, bool& success );

  bool hasOnlyCHneighbors(int atom ) const;
  bool hasAtLeastOneNeighborOfElement(int atom, const std::string& elem ) const;
  bool hasExactlyTwoOxygenNeighbors(int atom ) const;
  int countHneighbors(int atom) const;

};

#endif /*ATOMCATEGORIES_H_*/
