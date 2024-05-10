#include <cmath>
#include "EMStructure.h"


////////////////////////////////////////////////////////////////////////////////
// Description:
//   Returns the square of the distance between a point in Cartesian
//   space and the location of an atom.
////////////////////////////////////////////////////////////////////////////////
double EMStructure::distanceSquared(int id, Vec3 pos) {
  Vec3 temp = rigidUnitSystem->meanPositions(id);
  return temp.dist2(pos);
}


////////////////////////////////////////////////////////////////////////////////
// Description:
//   Calculates the unscaled density contributed by an atom at a distance r.  
//   This density will later be scaled by the mass of the contributing atom. 
// Parameters:
//   double r2 -- the square of the distance from the atom
// Return Value List:
//   Returns the calculated contribution to the electron density at r.
////////////////////////////////////////////////////////////////////////////////
double EMStructure::scaledDensity(double r2, int atomID) {
  if (r2 > cutoff*cutoff) { return 0.0; }
  int bin = (int) (r2 / binSize);
  return lookupTable[bin] * atoms[atomID].Z;
}

///////////////////////////////////////////////////////////////////////////////
// Description:
//   Constructor for the EMStructure class.  Creates an EMStructure object 
//   that will be associated with conformers produced by newsim.
// Parameters:
//   PDB &pdb -- an object of type PDB that was created when the starting
//           structure was reas in by newsim; used to get atomic charges
//   RigidUnitSystem &RUSInput -- This will contain (among other things) the
//           element names of the Atom objects contained in the 
//           EMStructure::atoms vector
//   double resolution -- resolution factor used by 
//           EMStructure::densityAt() 
//   double contained -- the fraction of the density contribution from each
//           atom that should be included in the cutoff.  The approximation
//           used below to get the cutoff is good for values > 0.3 and < 1.0
///////////////////////////////////////////////////////////////////////////////
EMStructure::EMStructure(PDB &pdb, RigidUnitSystem &RUSInput, 
                         double resolution, double contained) {
  sigma = 0.5 * resolution; 
  // NOTE: Situs uses 0.5*resolution, Topf et al. Structure 16:295-307 (2008)
  // uses 0.425*resolution.  
  double a = 0.249*(log(36) + log(contained/(1-contained)));
  cutoff = a*sigma;
  rigidUnitSystem = &RUSInput;
  grid = new Grid(*rigidUnitSystem,pdb,cutoff);

  // set up lookup table
  binSize = 0.01; //default value
  size_t numBins = (unsigned int) (cutoff*cutoff / binSize);
  lookupTable.resize(numBins+1);
  for (size_t i = 0; i <= numBins; i++) {
    double r2 = i * binSize;
    lookupTable[i] = exp(-1.5*r2/(sigma*sigma));
  }
  // set up atoms
  for (size_t p = 0; p < pdb.atomLines.size(); p++) {
    short int Z;
    if (pdb.atomLines[p].element == "C") { Z = 6; }
    else if (pdb.atomLines[p].element == "N") { Z = 7; }
    else if (pdb.atomLines[p].element == "O") { Z = 8; }
    else if (pdb.atomLines[p].element == "P") { Z = 15; }
    else if (pdb.atomLines[p].element == "S") { Z = 16; }
    else if (pdb.atomLines[p].element == "MG") { Z = 12; }
    else if (pdb.atomLines[p].element == "MN") { Z = 25; }
    else if (pdb.atomLines[p].element == "SI") { Z = 14; }
    else if (pdb.atomLines[p].element == "FE") { Z = 26; }
    else if (pdb.atomLines[p].element == "H") { Z = 0; }
    // hydrogens are added here, but not in the coarse grid, so that atom numbers here are the
    // same as in rigidUnitSystem
    else { Z = 6; } // assume carbon
    // TODO:  I should probably make sure that I can do this -- I'm assuming that the ordering
    // of points in the RigidUnitSystem is the same as the ordering of lines in pdb.atomLines[].
    Atom newAtom;
    newAtom.id = p;
    newAtom.Z = Z;
    atoms.push_back(newAtom);
  }
    
  // the old FRODA Grid class isn't here anymore; I'll probably want to introduce
  // some kind of coarse grid

  return;
}

///////////////////////////////////////////////////////////////////////////////
// Description:
//   Destructor for the EMStructure class. 
///////////////////////////////////////////////////////////////////////////////
EMStructure::~EMStructure() {
  delete grid;
}

///////////////////////////////////////////////////////////////////////////////
// Description:
//   Calculates the density at a given point in space resulting from 
//   the atoms contained in an EMStructure object.
// Parameters:
//   Vec3 v -- the location in Cartesian space where the density is to be
//               calculated
// Return Value List:
//   Returns the calculated density
///////////////////////////////////////////////////////////////////////////////
double EMStructure::densityAt(Vec3 v) {
  using namespace std;
  double totalDensity = 0.0;
  // this version doesn't use the coarse grid
  /*for (size_t i = 0; i < atoms.size(); i++) {
    double r2 = distanceSquared(atoms[i].id,v);
    if (r2 < cutoff*cutoff) {
      totalDensity += unscaledDensity(r2) * atoms[i].Z;    
    } 
  }*/
  std::vector<int> nearbyAtoms = grid->nearbyAtoms(v);
  for (size_t i = 0; i < nearbyAtoms.size(); i++) {
    double r2 = distanceSquared(atoms[nearbyAtoms[i]].id,v);
    if (r2 < cutoff*cutoff) {
      totalDensity += scaledDensity(r2,nearbyAtoms[i]);
    }  
  }
  //cerr << totalDensity << endl;
  return totalDensity;
}
