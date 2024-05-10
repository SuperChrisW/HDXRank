/*
 * TargetingConstraintManager.cpp
 *
 *  Created on: Sep 17, 2009
 *      Author: dwfarrel
 */

#include "ProteinInfo.h"
#include "RigidUnitSystem.h"
#include "ConstraintEnforcingPotential_Targeting.h"
#include "Settings.h"
#include "ObjectBuilder_FromSettings.h"
#include "TextFileInput.h"
#include <iomanip>
#include <sstream>
#include <algorithm>

using namespace std;

ConstraintEnforcingPotential_Targeting::ConstraintEnforcingPotential_Targeting(
  RigidUnitSystem* sys,
  const ProteinInfo& prot,
  const NeighborTable& nt,
  const Settings& settings ) :
    ConstraintEnforcingPotential( sys, prot, nt, settings ),
    targetEnergy( NULL ),
    targmap( NULL ),
    targneigh( NULL ),
    fixedCommon( false ),
    fixedNoncommon( false ),
    removeCommon( false ),
    removeNoncommon( false )
{
  if ( settings.outputconstraintlists ) {
    if ( bbhb ) bbhb->write("initial_bbhb.txt");
    if ( hb ) hb->write("initial_hb.txt");
    if ( ph ) ph->write("initial_ph.txt");
  }

  //First, attempt to get optional data regarding the target structure.
  //These optional target data are necessary for merging constraints
  //between the two sets. It the user provides the data, th
  //specifically the proteinInfo object, the neighbor table,
  //the rigid unit system, and the constraint enforcing potential.
  //We can only get this data if the user has supplied enough information
  //in the settings.
  ProteinInfo *targprot = NULL;
  NeighborTable* targnt = NULL;
  RigidUnitSystem* targsys = NULL;
  ConstraintEnforcingPotential* targcep = NULL;

  if ( settings.mergeConstraintsForTargeting &&
         settings.files.targpdb != "" &&
         settings.files.targfirstcov != "" &&
         settings.files.targfirstrc != "" &&
         settings.files.targatomtypes != "" ) {
    Settings altersettings;
    altersettings = settings;
    altersettings.files.pdb = settings.files.targpdb;
    altersettings.files.firstcov = settings.files.targfirstcov;
    altersettings.files.firstrc = settings.files.targfirstrc;
    altersettings.files.atomTypesAssignment = settings.files.targatomtypes;
    altersettings.doRedetermineConstraints = false;
    altersettings.files.overrideMinDistFile = "";

    cout << "\nLoading Target Structure Data" << endl;
    ObjectBuilder_FromSettings* objBuilder = new ObjectBuilder_FromSettings( altersettings );
    targprot = objBuilder->getProteinInfo();
    targnt = objBuilder->getCovalentNeighborTable();
    targsys = objBuilder->getRigidUnitSystem();
    targcep = new ConstraintEnforcingPotential( targsys, *targprot, *targnt, altersettings );

    if ( settings.outputconstraintlists ) {
      if ( targcep->bbhb ) targcep->bbhb->write("target_bbhb.txt");
      if ( targcep->hb ) targcep->hb->write("target_hb.txt");
      if ( targcep->ph ) targcep->ph->write("target_ph.txt");
    }

    delete targcep->overlapEnergy;
    targcep->overlapEnergy = NULL;

    delete objBuilder;
  }

  //initialize the target map and the target energy term
  if ( targsys ) {
    initialize( sys, prot, settings, targsys->meanPositions() );
  }
  else {
    ObjectBuilder_FromSettings* objBuilder = new ObjectBuilder_FromSettings( settings );
    vector<Vec3> targcoords;
    cout << "Loading target coordinates" << endl;
    objBuilder->getTargetCoords( targcoords );

    initialize( sys, prot, settings, targcoords );
    delete objBuilder;
  }

  cout << endl;

  //merge the minimum distance constraints of the two structures,
  //taking the lower of the two constraints
  if ( targcep ) {
    //merge side chain torsion minimum distance constraints
    if ( sideChainTorsion && targcep->sideChainTorsion ) {
      cout << "Merging side-chain torsion constraints" << endl;
      setupConstraintsSC( sideChainTorsion, targcep->sideChainTorsion );
    }
    //apply the min-dist constraint overrides from the target structure
    //to the
    if ( settings.doAutoOverrideMinDist ) {
      cout << "Lowering any minimum distance constraints that are violated in target structure" << endl;
      addMinDistOverridesFromTarg( targcep );
    }
  }

  //Merge the maximum-distance constraints of the two structures,
  //taking the higher of the two constraints
  if ( targcep ) {
    cout << "Merging hydrogen bond and hydrophobic constraints" << endl;

    fixedCommon = settings.commonConstraintHandling == "fixed";
    fixedNoncommon = settings.noncommonConstraintHandling == "fixed";
    removeCommon = settings.commonConstraintHandling == "remove";
    removeNoncommon = settings.noncommonConstraintHandling == "remove";

    cout << " Common Hydrogen Bond and Hydrophobic Contact Constraints: " << settings.commonConstraintHandling << endl;
    cout << " Noncommon Hydrogen Bond and Hydrophobic Contact Constraints: " << settings.noncommonConstraintHandling << endl;

    //merge hbonds
    if ( bbhb && targcep->bbhb ) {
      setupConstraintsBBHB( bbhb, targcep->bbhb );
    }
    if ( hb && targcep->hb ) {
      setupConstraintsHB( hb, targcep->hb );
    }
    //merge phobics
    if ( ph && targcep->ph ) {
      setupConstraintsPH( ph, targcep->ph );
    }
    notifyTermsChanged();
  }

  //We're done with the optional target data, so delete it all
  delete targcep;
  delete targsys;
  delete targprot;

  //save the target neighbor table
  if ( targnt ) {
    targneigh = new NeighborTable( targnt->size() );
    targneigh->insert( *targnt );
    targneigh->commit();
  }
  delete targnt;

}

ConstraintEnforcingPotential_Targeting::~ConstraintEnforcingPotential_Targeting() {
  delete targetEnergy;
  delete targneigh;
  delete targmap;
}

void ConstraintEnforcingPotential_Targeting::initialize(
  RigidUnitSystem* sys,
  const ProteinInfo& prot,
  const Settings& settings,
  const vector<Vec3>& targcoords ) {

  ObjectBuilder_FromSettings* objBuilder = new ObjectBuilder_FromSettings( settings );
  //setup target map
  targmap = objBuilder->getTargetMap();
  delete objBuilder;

  //if no atommap was provided, but the two states have the
  //same number of atoms, assume there is a trivial 1-to-1 mapping
  int natoms = sys->nPoints();
  if ( !targmap && targcoords.size() == static_cast<size_t>( natoms ) ) {
    targmap = new TargetMap;
    for ( int i = 0; i < natoms; i++ ) {
      targmap->insert( i, i );
    }
  }
  else if ( !targmap ) {
    cout << "Error:  Missing atommap" << endl;
    exit(0);
  }

  set<int> targsubset;
  bool useTargSubset = false;
  if ( settings.files.targsubset != "" ) {
    TextFileInput t;
    t.readAtomSet( settings.files.targsubset, targsubset );
    useTargSubset = true;
  }

  if ( settings.files.breakableconstraints_atomlist != "" ) {
    TextFileInput t;
    t.readAtomSet( settings.files.breakableconstraints_atomlist, breakable_atomlist );
  }

  //instantiate global target energy and incorporate into energy/gradient
  cout << "Establishing RMSD-to-target energy term" << endl;
  targetEnergy = new TargetEnergy( sys );
  map<int,int>::const_iterator end = targmap->src2targ().end();
  for ( map<int,int>::const_iterator it = targmap->src2targ().begin(); it != end; it++ ) {
    if ( useTargSubset && targsubset.find( it->first ) == targsubset.end() ) continue;
    string elem = prot.atom( it->first ).elem();
    if ( settings.targSubsetHeavy && elem == "H" ) continue;
    targetEnergy->addTargetPoint( it->first, targcoords[it->second] );
  }
  targetEnergy->update();
  targetEnergy->disable();
  addEnergyTerm( targetEnergy );
  addGradientTerm( targetEnergy );

  notifyTermsChanged();

}


string ConstraintEnforcingPotential_Targeting::generateColumnHeaderString() {
  ostringstream o;
  int i = 1;
  o << "Columns:\n";
  if ( sharedPointsEnergy ) o << i++ << " Worst Shared-Point Distance\n";
  if ( overlapEnergy || overridingMinDist )  o << i++ << " Worst Overlap Distance\n";
  if ( bbhb && hb ) o << i++ << " Worst H-bond distance violation\n";
  if ( bbhb && hb ) o << i++ << " Worst H-bond angle violation (deg)\n";
  if ( ph ) o << i++ << " Worst hydrophobic distance violation\n";
  if ( rama || overridingRama ) o << i++ << " Worst Rama distance violation\n";
  if ( sideChainTorsion || overridingSC ) o << i++ << " Worst side-chain torsion distance violation\n";
  o << i++ << " RMSD constraint\n";
  o << i++ << " RMSD to target\n";
  o << i++ << " projection along linear target coordinate\n";
  o << i++ << " projection along all orthogonal coordinates\n";
  return o.str();
}

string ConstraintEnforcingPotential_Targeting::generateSummaryString() {
  ostringstream o;
  o << fixed;
  if ( sharedPointsEnergy ) o << " " << setw(8) << sharedPointsEnergy->mismatch();
  if ( overlapEnergy || overridingMinDist ) o << " " << setw(8) <<
    max( overlapEnergy ? overlapEnergy->mismatch() : 0,
         overridingMinDist ? overridingMinDist->worstDistanceViolation() : 0 );
  if ( bbhb && hb ) o << " " << setw(8) << max( bbhb->worstDistanceViolation(), hb->worstDistanceViolation() );
  if ( bbhb && hb ) o << " " << setw(8) << max( bbhb->worstAngleViolationDeg(), hb->worstAngleViolationDeg() );
  if ( ph ) o << " " << setw(8) << ph->worstDistanceViolation();
  if ( rama || overridingRama ) o << " " << setw(8) <<
    max( rama ? rama->worstDistanceViolation() : 0,
         overridingRama ? overridingRama->worstDistanceViolation() : 0 );
  if ( sideChainTorsion || overridingSC ) o << " " << setw(8) <<
    max( sideChainTorsion ? sideChainTorsion->worstDistanceViolation() : 0,
         overridingSC ? overridingSC->worstDistanceViolation() : 0 );
  double coord1, coord2, rmsd;
  targetEnergy->getProjections( coord1, coord2, rmsd );
  //if ( code == 'A' ) o << " " << setw(10) << targetEnergy->getRMSDconstraint();
  //else               o << " " << setw(10) << "";
  o << " " << setw(10) << targetEnergy->getRMSDconstraint();
  o << " " << setw(10) << rmsd <<
          " " << setw(10) << coord1 <<
          " " << setw(10) << coord2;

  o << '\n';
  return o.str();
}

void ConstraintEnforcingPotential_Targeting::addMinDistOverridesFromTarg(
  const ConstraintEnforcingPotential* ceptarg  ) {

  int targp1;
  int targp2;
  int p1;
  int p2;
  bool success1;
  bool success2;
  MinDistConstraintContainer::iterator it;

  //whatever constraints need to be overridden in the target structure,
  //shall be overridden in the current structure.
  if ( ceptarg->overridingMinDist ) {
    for ( it = ceptarg->overridingMinDist->begin(); it != ceptarg->overridingMinDist->end(); it++ ) {
      targp1 = it->getp1();
      targp2 = it->getp2();
      targmap->conv2To1( targp1, p1, success1 );
      targmap->conv2To1( targp2, p2, success2 );
      if ( success1 && success2 ) overrideConstraint( p1, p2, it->getCutoff() );
    }
  }
  if ( ceptarg->overridingSC ) {
    for ( it = ceptarg->overridingSC->begin(); it != ceptarg->overridingSC->end(); it++ ) {
      targp1 = it->getp1();
      targp2 = it->getp2();
      targmap->conv2To1( targp1, p1, success1 );
      targmap->conv2To1( targp2, p2, success2 );
      if ( success1 && success2 ) overrideConstraint( p1, p2, it->getCutoff() );
    }
  }
  if ( ceptarg->overridingRama ) {
    for ( it = ceptarg->overridingRama->begin(); it != ceptarg->overridingRama->end(); it++ ) {
      targp1 = it->getp1();
      targp2 = it->getp2();
      targmap->conv2To1( targp1, p1, success1 );
      targmap->conv2To1( targp2, p2, success2 );
      if ( success1 && success2 ) overrideConstraint( p1, p2, it->getCutoff() );
    }
  }

  if ( overlapEnergy ) overlapEnergy->refresh();
  notifyTermsChanged();
}

void ConstraintEnforcingPotential_Targeting::setupConstraintsBBHB(
  BBHBContainer* bbhb, BBHBContainer* targbbhb ) {

  int h;
  int a;
  int targh;
  int targa;
  bool successh;
  bool successa;

  int nInitial = bbhb->size();
  int nCommon = 0;

  //For each constraint in structure 1, check to see if it is also found in structure 2.
  //If not found, then leave it as breakable.
  //If found, then set it as unbreakable.  Also inspect the distances
  //and angles and make them the most lenient of the two.
  BBHBContainer::iterator it = bbhb->begin();
  while ( it != bbhb->end() ) {
    //get the h and a for this constraint
    h = it->geth();
    a = it->geta();

    //lookup the h,a pair in the atommap to find the corresponding atom indices
    //in the target structure
    targmap->conv1To2( h, targh, successh );
    targmap->conv1To2( a, targa, successa );

    BBHBContainer::iterator ittarg;
    bool common;
    if ( successh && successa ) {
      ittarg = targbbhb->find( targh, targa );
      common = ittarg != targbbhb->end();
    }
    else {
      ittarg = targbbhb->end();
      common = false;
    }

    if ( common ) {
      if ( ittarg->getConstraintMaxDistHA() > it->getConstraintMaxDistHA() ) {
        it->setMaxDistHA( ittarg->getConstraintMaxDistHA() );
      }
      if ( ittarg->getConstraintMinAngleRadDHA() < it->getConstraintMinAngleRadDHA() ) {
        it->setMinAngleRadDHA( ittarg->getConstraintMinAngleRadDHA() );
      }
      if ( ittarg->getConstraintMinAngleRadHAB() < it->getConstraintMinAngleRadHAB() ) {
        it->setMinAngleRadHAB( ittarg->getConstraintMinAngleRadHAB() );
      }
      nCommon++;
    }

    if ( common && fixedCommon || !common && fixedNoncommon ) {
      if ( breakable_atomlist.find( h ) == breakable_atomlist.end() &&
        breakable_atomlist.find( a ) == breakable_atomlist.end() ) {
        it->makeUnbreakable();
      }
    }

    if ( common && removeCommon || !common && removeNoncommon ) {
      bbhb->erase( it );
      //DO NOT advance the iterator "it".  The erase operation
      //of ConstraintContainer replaces the current constraint
      //with the final constraint.  So, the "next" constraint
      //to be checked is the one that "it" is now currently pointing to.
    }
    else it++;
  }

  cout << "Backbone-backbone Hydrogen Bond Constraints: " << endl;
  cout << "  Initial State: " << nInitial << endl;
  cout << "  Target State: " << targbbhb->size() << endl;
  cout << "  Common: " << nCommon << endl;
  cout << "  Initial State Non-common: " << nInitial - nCommon << endl;

}

void ConstraintEnforcingPotential_Targeting::setupConstraintsHB(
  HBContainer* hb, HBContainer* targhb ) {
  int h;
  int a;
  int targh;
  int targa;
  bool successh;
  bool successa;

  int nInitial = hb->size();
  int nCommon = 0;

  HBContainer::iterator it = hb->begin();
  while ( it != hb->end() ) {
    //get the h and a for this constraint
    h = it->geth();
    a = it->geta();

    //lookup the h,a pair in the atommap to find the corresponding atom indices
    //in the target structure
    targmap->conv1To2( h, targh, successh );
    targmap->conv1To2( a, targa, successa );

    HBContainer::iterator ittarg;
    bool common;
    if ( successh && successa ) {
      ittarg = targhb->find( targh, targa );
      common = ittarg != targhb->end();
    }
    else {
      ittarg = targhb->end();
      common = false;
    }

    if ( common ) {
      if ( ittarg->getConstraintMaxDistHA() > it->getConstraintMaxDistHA() ) {
        it->setMaxDistHA( ittarg->getConstraintMaxDistHA() );
      }
      if ( ittarg->getConstraintMinAngleRadDHA() < it->getConstraintMinAngleRadDHA() ) {
        it->setMinAngleRadDHA( ittarg->getConstraintMinAngleRadDHA() );
      }
      nCommon++;
    }

    if ( common && fixedCommon || !common && fixedNoncommon ) {
      if ( breakable_atomlist.find( h ) == breakable_atomlist.end() &&
        breakable_atomlist.find( a ) == breakable_atomlist.end() ) {
        it->makeUnbreakable();
      }
    }

    if ( common && removeCommon || !common && removeNoncommon ) {
      hb->erase( it );
      //DO NOT advance the iterator "it".  The erase operation
      //of ConstraintContainer replaces the current constraint
      //with the final constraint.  So, the "next" constraint
      //to be checked is the one that "it" is now currently pointing to.
    }
    else it++;
  }

  cout << "All other Hydrogen Bond Constraints: " << endl;
  cout << "  Initial State: " << nInitial << endl;
  cout << "  Target State: " << targhb->size() << endl;
  cout << "  Common: " << nCommon << endl;
  cout << "  Initial State Non-common: " << nInitial - nCommon << endl;

}

void ConstraintEnforcingPotential_Targeting::setupConstraintsPH(
  PHContainer* ph, PHContainer* targph ) {
  int p1;
  int p2;
  int targp1;
  int targp2;
  bool success1;
  bool success2;

  int nInitial = ph->size();
  int nCommon = 0;

  PHContainer::iterator it = ph->begin();
  while ( it != ph->end() ) {
    //get the two atoms for this constraint.
    p1 = it->getp1();
    p2 = it->getp2();

    //lookup the p1,p2 pair in the atommap to find the corresponding atom indices
    //in the target structure
    targmap->conv1To2( p1, targp1, success1 );
    targmap->conv1To2( p2, targp2, success2 );

    PHContainer::iterator ittarg;
    bool common;
    if ( success1 && success2 ) {
      ittarg = targph->find( targp1, targp2 );
      common = ittarg != targph->end();
    }
    else {
      ittarg = targph->end();
      common = false;
    }

    if ( common ) {
      if ( ittarg->getCutoff() > it->getCutoff() ) {
        it->setCutoff( ittarg->getCutoff() );
      }
      nCommon++;
    }

    if ( common && fixedCommon || !common && fixedNoncommon ) {
      if ( breakable_atomlist.find( p1 ) == breakable_atomlist.end() &&
        breakable_atomlist.find( p2 ) == breakable_atomlist.end() ) {
        it->makeUnbreakable();
      }
    }

    if ( common && removeCommon || !common && removeNoncommon ) {
      ph->erase( it );
      //DO NOT advance the iterator "it".  The erase operation
      //of ConstraintContainer replaces the current constraint
      //with the final constraint.  So, the "next" constraint
      //to be checked is the one that "it" is now currently pointing to.
    }
    else it++;
  }

  cout << "Hydrophobic Constraints: " << endl;
  cout << "  Initial State: " << nInitial << endl;
  cout << "  Target State: " << targph->size() << endl;
  cout << "  Common: " << nCommon << endl;
  cout << "  Initial State Non-common: " << nInitial - nCommon << endl;

}

void ConstraintEnforcingPotential_Targeting::setupConstraintsSC(
  SideChainTorsionContainer* tor,
  SideChainTorsionContainer* targtor ) {

  int p1;
  int p2;
  int targp1;
  int targp2;
  bool success1;
  bool success2;

  for ( SideChainTorsionContainer::iterator it = tor->begin();
        it != tor->end(); it++ ) {
    //get the two atoms for this constraint.
    p1 = it->getp1();
    p2 = it->getp2();

    //lookup the p1,p2 pair in the atommap to find the corresponding atom indices
    //in the target structure
    targmap->conv1To2( p1, targp1, success1 );
    targmap->conv1To2( p2, targp2, success2 );
    if ( success1 && success2 ) {
      //now see if the constraint is held in common.
      SideChainTorsionContainer::iterator ittarg;
      ittarg = targtor->find( targp1, targp2 );
      if ( ittarg != targtor->end() ) {
        //use whichever cutoff is SMALLER in this case
        if ( ittarg->getCutoff() < it->getCutoff() ) {
          it->setCutoff( ittarg->getCutoff() );
        }
      }
    }
  }
}
