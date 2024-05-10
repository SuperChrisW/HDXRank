#include "ObjectBuilder_FromSettings.h"
#include "RigidUnitSystemBuilder.h"
#include "RigidUnitSystem.h"
#include "ConstraintEnforcingPotential.h"
#include "PerturbRelaxCycle.h"
#include "Settings.h"
#include "Observable.h"
#include "ProteinInfo.h"
#include "PDB.h"
#include "AmberPrmtop.h"
#include "TextFileInput.h"
#include "FIRSTFileInput.h"
#include "AmberRestartFileInput.h"
#include "TargetMap.h"
#include <string>
#include <iostream>

using namespace std;

ObjectBuilder_FromSettings::ObjectBuilder_FromSettings(
  const Settings &settings_ ) :
    settings(settings_),
    pdb( NULL ),
    prmtop( NULL ),
    targpdb( NULL ),
    targprmtop( NULL )
{
}

ObjectBuilder_FromSettings::~ObjectBuilder_FromSettings()
{
  delete pdb;
  delete prmtop;
  delete targpdb;
  delete targprmtop;
}

void ObjectBuilder_FromSettings::initializeInputData() {
  if ( !pdb && settings.files.pdb != "" ) pdb = new PDB( settings.files.pdb );
  if ( !prmtop && settings.files.prmtop != "" ) prmtop = new AmberPrmtop( settings.files.prmtop );
}

void ObjectBuilder_FromSettings::initializeTargetData() {
  if ( !targpdb && settings.files.targpdb != "" ) targpdb = new PDB( settings.files.targpdb );
  if ( !targprmtop && settings.files.targprmtop != "" ) targprmtop = new AmberPrmtop( settings.files.targprmtop );
}

ProteinInfo* ObjectBuilder_FromSettings::getProteinInfo() {
  initializeInputData();

  if ( pdb ) return new ProteinInfo( *pdb );
  else if ( prmtop ) return new ProteinInfo( *prmtop );
  else {
    cout << "Error: Could not build ProteinInfo object \n";
    cout << "       because pdb or prmtop was not specified." << endl;
    exit(0);
  }
  return NULL;
}
ProteinInfo* ObjectBuilder_FromSettings::getTargetProteinInfo() {
  initializeTargetData();
  if ( targpdb ) return new ProteinInfo( *targpdb );
  else if ( targprmtop ) return new ProteinInfo( *targprmtop );
  else {
    cout << "Error: Could not build ProteinInfo object \n";
    cout << "       because targpdb or targprmtop was not specified." << endl;
    exit(0);
  }
  return NULL;
}

NeighborTable* ObjectBuilder_FromSettings::getCovalentNeighborTable() {
  initializeInputData();
  NeighborTable* nt = NULL;
  if ( prmtop ) {
    TextFileInput textFileInput;
    nt = textFileInput.buildCovalentNeighborTable_FromPrmtop( *prmtop );
  }
  else if ( pdb && settings.files.firstcov != "" ) {
    FIRSTFileInput firstFileInput( pdb, settings.files.FIRSTBondFileInterpretation );
    nt = firstFileInput.buildCovalentNeighborTable_FromFIRSTcov( settings.files.firstcov );
  }
  else {
    cout << "Error: covalent neighbor table not specified" << endl;
    exit(0);
  }
  return nt;
}

NeighborTable* ObjectBuilder_FromSettings::getTargetCovalentNeighborTable() {
  initializeTargetData();
  NeighborTable* nt = NULL;
  if ( targprmtop ) {
    TextFileInput textFileInput;
    nt = textFileInput.buildCovalentNeighborTable_FromPrmtop( *targprmtop );
  }
  else if ( targpdb && settings.files.targfirstcov != "" ) {
    FIRSTFileInput firstFileInput( targpdb, settings.files.FIRSTBondFileInterpretation );
    nt = firstFileInput.buildCovalentNeighborTable_FromFIRSTcov( settings.files.targfirstcov );
  }
  return nt;
}

TargetMap* ObjectBuilder_FromSettings::getTargetMap() {
  TargetMap* atommap = NULL;
  TextFileInput textFileInput;
  if ( settings.files.atommap != "") {
    atommap = textFileInput.buildTargetMap( settings.files.atommap );
  }
  return atommap;

}

void ObjectBuilder_FromSettings::getTargetCoords( vector<Vec3>& targCoords ) {
  initializeTargetData();
  if ( targpdb ) targpdb->getPositions( targCoords );
  else if ( targprmtop && settings.files.targAmberRestart != "" ) {
    AmberRestartFileInput amb;
    bool good;
    amb.read( settings.files.targAmberRestart, *targprmtop, targCoords, good );
    if ( !good ) targCoords.clear();
  }
}

RigidUnitSystem* ObjectBuilder_FromSettings::getRigidUnitSystem() {
  initializeInputData();
  RigidUnitSystemBuilder *builder = new RigidUnitSystemBuilder();
  builder->setpdbfilename( settings.files.pdb );
  builder->setprmtopfilename( settings.files.prmtop );
  builder->setfirstrcfilename( settings.files.firstrc );
  builder->setfirstcovfilename( settings.files.firstcov );
  builder->setfirsthbondsfilename( settings.files.firsthbonds, settings.files.hbondEnergyCutoff );
  builder->setFIRSTBondFileInterpretation( settings.files.FIRSTBondFileInterpretation );
  builder->setrestartpdbfilename( settings.files.restartpdb );

  RigidUnitSystem* sys = builder->build();

  delete builder;

  return sys;
}
