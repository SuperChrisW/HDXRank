#include "AtomCategories.h"
#include "ProteinInfo.h"
#include "NeighborTable.h"
#include <iostream>

using namespace std;

AtomCategories::AtomCategories( const NeighborTable *neighborTable_, const ProteinInfo *prot_ ) :
  neighborTable( neighborTable_ ),
  prot( prot_ )
{
  defineTypes();

  //only if residue matches one of these will the
  //atom be assigned to a protein atom type
  proteinResidueNames.insert( "LEU" );
  proteinResidueNames.insert( "ALA" );
  proteinResidueNames.insert( "VAL" );
  proteinResidueNames.insert( "MET" );
  proteinResidueNames.insert( "ILE" );
  proteinResidueNames.insert( "TRP" );
  proteinResidueNames.insert( "PHE" );
  proteinResidueNames.insert( "SER" );
  proteinResidueNames.insert( "THR" );
  proteinResidueNames.insert( "CYS" );
  proteinResidueNames.insert( "CYX" ); //amber name
  proteinResidueNames.insert( "ASN" );
  proteinResidueNames.insert( "GLN" );
  proteinResidueNames.insert( "TYR" );
  proteinResidueNames.insert( "HIS" );
  proteinResidueNames.insert( "HID" ); //amber name
  proteinResidueNames.insert( "HIE" ); //amber name
  proteinResidueNames.insert( "ASP" );
  proteinResidueNames.insert( "GLU" );
  proteinResidueNames.insert( "LYS" );
  proteinResidueNames.insert( "ARG" );
  proteinResidueNames.insert( "GLY" );
  proteinResidueNames.insert( "PRO" );

  aromaticCarbonNames.insert("ARGCZ");
  aromaticCarbonNames.insert("PHECG");
  aromaticCarbonNames.insert("PHECD1");
  aromaticCarbonNames.insert("PHECE1");
  aromaticCarbonNames.insert("PHECZ");
  aromaticCarbonNames.insert("PHECE2");
  aromaticCarbonNames.insert("PHECD2");
  aromaticCarbonNames.insert("TRPCZ2");
  aromaticCarbonNames.insert("TRPCH2");
  aromaticCarbonNames.insert("TRPCZ3");
  aromaticCarbonNames.insert("TRPCE3");
  aromaticCarbonNames.insert("TYRCG");
  aromaticCarbonNames.insert("TYRCD1");
  aromaticCarbonNames.insert("TYRCE1");
  aromaticCarbonNames.insert("TYRCE2");
  aromaticCarbonNames.insert("TYRCD2");

  amideNitrogenNames.insert("GLNNE2");
  amideNitrogenNames.insert("ASNND2");

  carboxylOxygenNames.insert("ASPOD1");
  carboxylOxygenNames.insert("ASPOD2");
  carboxylOxygenNames.insert("GLUOE1");
  carboxylOxygenNames.insert("GLUOE2");

  assignTypesToAtoms();
}

AtomCategories::~AtomCategories()
{
}

void AtomCategories::defineTypes() {
  //establish PROTEIN atom types
  types.push_back("protCTa"); //0
  types.push_back("protCTb");  //1
  types.push_back("protC");                 //2
  types.push_back("protCA");
  types.push_back("protCother");
  types.push_back("protN");                 //5
  types.push_back("protN3");
  types.push_back("protNother");
  types.push_back("protO");
  types.push_back("protO2");
  types.push_back("protOH");                //10
  types.push_back("protOother");
  types.push_back("protS");
  types.push_back("protHN");
  types.push_back("protHS");
  types.push_back("protHO");                //15
  types.push_back("protHA");
  types.push_back("protHC");
  types.push_back("protHother");            //18

  //non-protein atom types
  types.push_back("C");                 //19
  types.push_back("H");                 //20
  types.push_back("N");
  types.push_back("O");
  types.push_back("S");
  types.push_back("P");
  types.push_back("MG");
  types.push_back("MN");
  types.push_back("SI");
  types.push_back("FE");
  types.push_back("XX");                //29

  for ( size_t i = 0; i < types.size(); i++ ) {
    mapTypeNameToTypeIndex[ types[i] ] = i;
  }

}

bool AtomCategories::hasOnlyCHneighbors(int atom) const {
  for ( size_t i = 0; i < (*neighborTable)[atom].size(); i++ ) {
    int neigh = (*neighborTable)[atom][i];
    if ( prot->atom(neigh).elem() != "C" &&
         prot->atom(neigh).elem() != "H" ) return false;
  }
  return true;
}

int AtomCategories::countHneighbors(int atom) const {
  int count = 0;
  for ( size_t i = 0; i < (*neighborTable)[atom].size(); i++ ) {
    int neigh = (*neighborTable)[atom][i];
    if ( prot->atom(neigh).elem() == "H" ) count++;
  }
  return count;
}

bool AtomCategories::hasAtLeastOneNeighborOfElement(int atom, const string& elem ) const {
  for ( size_t i = 0; i < (*neighborTable)[atom].size(); i++ ) {
    int neigh = (*neighborTable)[atom][i];
    if ( prot->atom(neigh).elem() == elem ) return true;
  }
  return false;
}

bool AtomCategories::hasExactlyTwoOxygenNeighbors(int atom ) const {
  int count = 0;
  for ( size_t i = 0; i < (*neighborTable)[atom].size(); i++ ) {
    int neigh = (*neighborTable)[atom][i];
    if ( prot->atom(neigh).elem() == "O" ) count++;
  }
  return count == 2;
}

void AtomCategories::getProteinAtomType( int atom, int& atomtype, bool& success ) {

  string elem = prot->atom(atom).elem();
  string resname = prot->atom(atom).resi().name();
  string atomname = prot->atom(atom).name();
  string resPlusName = resname + atomname;

  const vector<int> *neighbors = &(*neighborTable)[atom];
  int nNeighbors = neighbors->size();

  string atomtypestring;

  //Carbons
  if ( elem == "C" ) {
    if ( nNeighbors == 3 && hasAtLeastOneNeighborOfElement( atom, "O" ) )
      atomtypestring = "protC";
    else if ( nNeighbors == 4 ) {
      if ( hasOnlyCHneighbors( atom ) ) {
        atomtypestring = "protCTa";
      }
      else {
        atomtypestring = "protCTb";
      }
    }
    else if ( nNeighbors == 3 && find(
          aromaticCarbonNames.begin(),
          aromaticCarbonNames.end(),
          resPlusName ) != aromaticCarbonNames.end() )
    {
      atomtypestring = "protCA";
    }
    else atomtypestring = "protCother";
  }

  //Hydrogens
  else if ( elem == "H" ) {
    if ( nNeighbors == 1 ) {
      int neigh = (*neighbors)[0];
      string neighborElem = prot->atom(neigh).elem();
      if ( neighborElem == "N" ) {
        atomtypestring = "protHN";
      }
      else if ( neighborElem == "S" ) {
        atomtypestring = "protHS";
      }
      else if ( neighborElem == "O" ) {
        atomtypestring = "protHO";
      }
      else if ( neighborElem == "C" ) {
        string neighborResPlusName = prot->atom(neigh).resi().name() + prot->atom(neigh).name();
        if ( find(
          aromaticCarbonNames.begin(),
          aromaticCarbonNames.end(),
          neighborResPlusName ) != aromaticCarbonNames.end() )
        {
          atomtypestring = "protHA";
        }
        else if ( hasOnlyCHneighbors( neigh ) ) {
          atomtypestring = "protHC";
        }
        else atomtypestring = "protHother";
      }
      else atomtypestring = "protHother";
    }
    else atomtypestring = "protHother";
  }

  //Nitrogen
  else if ( elem == "N" ) {
    if ( nNeighbors == 4 )
      atomtypestring = "protN3";
    else if ( atomname == "N" ||
        find( amideNitrogenNames.begin(),
              amideNitrogenNames.end(),
              resPlusName ) != amideNitrogenNames.end() )
    {
      atomtypestring = "protN";
    }
    else
      atomtypestring = "protNother";
  }

  //Oxygen
  else if ( elem == "O" ) {
    if ( hasAtLeastOneNeighborOfElement( atom, "H" ) )
      atomtypestring = "protOH";
    else if ( nNeighbors == 1 ) {
      int neigh = (*neighbors)[0];
      string neighborElem = prot->atom(neigh).elem();
      if ( neighborElem == "C" &&
          (*neighborTable)[neigh].size() == 3 &&
          hasExactlyTwoOxygenNeighbors( neigh ) ) {
        atomtypestring = "protO2";
      }
      else if ( neighborElem == "C" && hasAtLeastOneNeighborOfElement( neigh, "N" ) ) {
        atomtypestring = "protO";
      }
      else atomtypestring = "protOother";
    }
    else atomtypestring = "protOother";
  }

  //Sulfur
  else if ( elem == "S" )
    atomtypestring = "protS";

  if ( atomtypestring == "" ) {
    atomtype = -1;
    success = false;
    cout << "Warning, could not find a standard protein atom type for protein atom" << endl;
    cout << resname << " " << atomname << endl;
  }
  else {
    map<string,int>::iterator it = mapTypeNameToTypeIndex.find( atomtypestring );
    if ( it != mapTypeNameToTypeIndex.end() ) {
      atomtype = it->second;
      success = true;
    }
    else {
      atomtype = -1;
      success = false;
      cout << "Error, atom type " << atomtypestring << " has not been assigned a type number" << endl;
      exit(0);
    }
  }

}

void AtomCategories::assignTypesToAtoms() {
  int natoms = prot->natoms();
  atomIndexToType.resize( natoms );

  int countProteinAtoms = 0;
  int countNonproteinAtoms = 0;
  int countUnrecognizedAtoms = 0;
  for ( int atom = 0; atom < natoms; atom++ ) {
    int atomtype = -1;
    bool success = false;

    //protein atoms
    if ( proteinResidueNames.find( prot->atom(atom).resi().name() ) != proteinResidueNames.end() ) {
      getProteinAtomType( atom, atomtype, success );
      if ( success ) countProteinAtoms++;
    }

    //non-protein atoms
    if ( !success ) {
      map<string, int>::const_iterator it = mapTypeNameToTypeIndex.find( prot->atom(atom).elem() );
      if ( it != mapTypeNameToTypeIndex.end() ) {
        atomtype = it->second;
        success = true;
        countNonproteinAtoms++;
      }
    }

    //Unrecognized type
    if ( !success ) {
      map<string,int>::iterator it = mapTypeNameToTypeIndex.find( "XX" );
      if ( it != mapTypeNameToTypeIndex.end() ) {
        atomtype = it->second;
        cout << "Warning: unrecognized atom " << atom << ", elem " << prot->atom(atom).elem()  <<
                ", type set to Unknown (type " << atomtype << ")" << endl;
        success = true;
        countUnrecognizedAtoms++;
      }
    }

    if ( !success ) {
      cout << "Error in atom type assignment" << endl;
      exit(0);
    }

    atomIndexToType[atom] = atomtype;
  }

  cout << countProteinAtoms << " protein atom types" << endl;
  cout << countNonproteinAtoms << " non-protein atom types" << endl;
  cout << countUnrecognizedAtoms << " unrecognized atom types" << endl;
}

