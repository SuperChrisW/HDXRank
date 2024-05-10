#ifndef GENERATEOVERLAPPINGRIGIDUNITS_H_
#define GENERATEOVERLAPPINGRIGIDUNITS_H_

class NeighborTable;
#include <vector>

void generateOverlappingRigidUnits( 
    const std::vector<int>& mapPtoRC,
    const NeighborTable& neighborTable,
    std::vector< std::vector<int> >& mapRUtoP_vector);

#endif /*GENERATEOVERLAPPINGRIGIDUNITS_H_*/
