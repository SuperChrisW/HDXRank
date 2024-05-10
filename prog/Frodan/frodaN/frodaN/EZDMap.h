// Author name:  Craig Jolley
// Created:      19 Feb 2008

// NOTE: Back in the days of classic FRODA, this class was contained in GenericMap.h.

#ifndef EZD_MAP_H_
#define EZD_MAP_H_

#include <string>
#include "GenericMap.h"

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Derived class of GenericMap -- this is designed to hold data taken from an
//   electron density file in the EZD (E-Z Density) format.  It inherits the
//   general functionalities from GenericMap and adds a few that are specific
//   to this data format.
////////////////////////////////////////////////////////////////////////////////
class EZDMap : public GenericMap {
private:
  char comment[80];
  double scale;        // scaling factor
  gridPoint num;       // number of grid points
public:
  EZDMap(std::string fileName);   // constructor for EZD_map
  ~EZDMap() {             // destructor    
    delete [] gridData; }             
  void displayHeader();    // display file header info
};

#endif
