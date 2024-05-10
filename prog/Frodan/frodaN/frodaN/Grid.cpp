#include "Grid.h"
#include <vector>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
//  A grid system that keeps track of atoms in blocks of real space.  Like
//  the Grid class in FIRST, only far more awesome.
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description: 
//   Constructor for Grid class
////////////////////////////////////////////////////////////////////////////////
Grid::Grid(RigidUnitSystem &rigidUnitSystem_, PDB &pdb, double gridLength_) {
  rigidUnitSystem = &rigidUnitSystem_;
  gridLength = gridLength_;
  contributingAtoms.clear();
  for (size_t p = 0; p < pdb.atomLines.size(); p++) {
    if (pdb.atomLines[p].element != "H") {contributingAtoms.push_back(p);}
  }
  setupGrid();
  updateGrid();
  return;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//   Finds the minimum and maximum atom coordinates in rigidUnitSystem
////////////////////////////////////////////////////////////////////////////////
void Grid::findMinMax() {
  int np = rigidUnitSystem->nPoints();
  min = max = rigidUnitSystem->meanPositions(0);
  for (int p = 1; p < np; p++) {
    Vec3 v = rigidUnitSystem->meanPositions(p);
    if (v.x > max.x) { max.x = v.x; }
    else if (v.x < min.x) { min.x = v.x; }
    if (v.y > max.y) { max.y = v.y; }
    else if (v.y < min.y) { min.y = v.y; }
    if (v.z > max.z) { max.z = v.z; }
    else if (v.z < min.z) { min.z = v.z; }
  }
  return;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//   Sets up an empty myGrid with the appropriate dimensions, and loads 
//   adjacentPoints with information to simplify subsequent searches
////////////////////////////////////////////////////////////////////////////////
void Grid::setupGrid() {
  findMinMax();
  Vec3 pad(30,30,30); // pad to allow for atom movement
  min -= pad;   
  max += pad;
  Vec3 gridSize = max - min;
  sizeX = (int) ceil(gridSize.x / gridLength) + 2;
  sizeY = (int) ceil(gridSize.y / gridLength) + 2;
  sizeZ = (int) ceil(gridSize.z / gridLength) + 2;
  myGrid.resize(sizeX*sizeY*sizeZ);
  // now set up adjacentPoints array
  adjacentPoints.resize(sizeX*sizeY*sizeZ);
  for (int x = 0; x < sizeX; x++) {
    int xi, xf;
    if (x > 0) { xi = x-1; } else { xi = x; } 
    if (x < sizeX-1) { xf = x+1; } else { xf = x; }
    for (int y = 0; y < sizeY; y++) {
      int yi, yf;
      if (y > 0) { yi = y-1; } else { yi = y; }
      if (y < sizeY-1) { yf = y+1; } else { yf = y; }
      for (int z = 0; z < sizeZ; z++) {
        int zi, zf;
        if (z > 0) { zi = z-1; } else { zi = z; }
        if (z < sizeZ-1) { zf = z+1; } else { zf = z; }
        int n = x + y*sizeX + z*sizeX*sizeY;
        for (int x2 = xi; x2 <= xf; x2++) {
          for (int y2 = yi; y2 <= yf; y2++) {
            for (int z2 = zi; z2 <= zf; z2++) {
              adjacentPoints[n].push_back(x2 + y2*sizeX + z2*sizeX*sizeY);
	    }
	  }
	}
      }
    }
  }
  return;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//   For a given position in real space, returns the 1-D array index of the 
//   grid point to which it belongs.
//////////////////////////////////////////////////////////////////////////////// 
int Grid::VecToArrayIndex(Vec3 v) {
  v -= min;
  int x = (int) (v.x / gridLength);
  int y = (int) (v.y / gridLength);
  int z = (int) (v.z / gridLength);
  // correct for out-of-bounds values
  if (x < 0) { x = 0; }
  if (x >= sizeX) { x = sizeX - 1; }
  if (y < 0) { y = 0; }
  if (y >= sizeY) { y = sizeY - 1; }
  if (z < 0) { z = 0; }
  if (z >= sizeZ) { z = sizeZ - 1; }
  return x + y*sizeX + z*sizeX*sizeY;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//   Updates grid contents based on atom positions in rigidUnitSystem
////////////////////////////////////////////////////////////////////////////////  
void Grid::updateGrid() {
  bool gridReset = false;
  do {
    // clear current contents of grid
    std::vector<std::vector<int> >::iterator igrid, endgrid;
    endgrid = myGrid.end();
    igrid = myGrid.begin();
    while (igrid != endgrid) {
      (*igrid++).clear();
    }
    // place atom numbers in grid
    gridReset = false;
    for (size_t i = 0; i < contributingAtoms.size(); i++) {
      int p = contributingAtoms[i];
      Vec3 v = rigidUnitSystem->meanPositions(p);
      if (v.x > max.x || v.y > max.y || v.z > max.z ||
          v.x < min.x || v.y < min.y || v.z < min.z) {
	std::cout << "WARNING: Atom is outside the EM coarse grid.  Re-configuring grid.\n";
        gridReset = true;
        break;
      }
      myGrid[VecToArrayIndex(v)].push_back(p);
    }
    if (gridReset) { setupGrid(); }
  } while (gridReset);
  return;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//   Returns atoms that are either in the same grid point as v, or in adjacent
//   ones.  This guarantees that all atoms that are within gridLength of v
//   will be in the list.
////////////////////////////////////////////////////////////////////////////////  
std::vector <int> Grid::nearbyAtoms(Vec3 v) {
  using namespace std;
  int n = VecToArrayIndex(v);
  vector <int> result;
  for (size_t i = 0; i < adjacentPoints[n].size(); i++) {
    int neighbor = adjacentPoints[n][i];
    result.insert(result.end(),myGrid[neighbor].begin(),myGrid[neighbor].end());
  }
  return result;
}

  
