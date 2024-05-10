#ifndef RIGIDUNITPOINTPDBBUILDER_H_
#define RIGIDUNITPOINTPDBBUILDER_H_

class PDB;
class RigidUnitSystem;
class NeighborTable;

class RigidUnitPointPDBBuilder
{
public:
  RigidUnitPointPDBBuilder(const PDB *pdb, const RigidUnitSystem *sys );
  virtual ~RigidUnitPointPDBBuilder();
  PDB *getPDB() { return newpdb; }
  //void addCONECTrecords( const NeighborTable& table );
private:
  PDB *newpdb;
};

#endif /*RIGIDUNITPOINTPDBBUILDER_H_*/
