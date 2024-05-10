/*
 * HBManager.cpp
 *
 *  Created on: Jun 11, 2009
 *      Author: dwfarrel
 */

#include "HBManager.h"
#include "FIRSTHbond.h"
#include "ProteinInfo.h"
#include "ForbidList.h"
#include <cmath>
#include <limits>
#include "mt19937ar.h"

const double deg2rad = M_PI/180.0;
const double rad2deg = 180.0/M_PI;

HBManager::HBManager(
  const ProteinInfo& proteinInfo,
  const NeighborTable& nt_,
  const vector<Vec3>& coords_,
  BBHBContainer& bb_hbonds_,
  HBContainer& hbonds_ ) :
    nt(nt_),
    coords(coords_),
    bb_hbonds( bb_hbonds_ ),
    hbonds( hbonds_ ),
    k( 10.0 ),
    L(2.0),
    theta2(140*deg2rad),
    theta1(100*deg2rad),
    nAdded_(0),
    nRemoved_(0),
    doAngleConstraints( false ),
    Ecutoff( -1.0 ),
    forbidlist( NULL ),
    carefulChecking( true )
{
  firstHbond = new FIRSTHbond( proteinInfo, nt, coords );

  //classify atoms
  int natoms = coords.size();
  isDonorH.assign( natoms, 0 );
  isAcceptor.assign( natoms, 0 );
  isBB.assign( natoms, 0 );
  for ( int i = 0; i < natoms; i++ ) {
    if ( firstHbond->isDonorH( i ) )
      isDonorH[i] = 1;
    else if ( firstHbond->isAcceptor( i ) )
      isAcceptor[i] = 1;

    if ( proteinInfo.atom(i).elem() == "H" && nt[i].size() == 1 && proteinInfo.atom( nt[i][0] ).name() == "N" ||
         proteinInfo.atom(i).name() == "N" ||
         proteinInfo.atom(i).name() == "O" ) {
      isBB[i] = 1;
    }
  }

  //prepare the verlet list
  verletCutoffs = new VerletCutoffsHAPairs( isDonorH, isAcceptor, &nt );
  verletCollector = new VerletCollectorHAPairs( isDonorH, isAcceptor );
  verlet = new VerletList( verletCutoffs, verletCollector, 1.0 );
  for ( int i = 0; i < natoms; i++ ) {
    if ( isDonorH[i] || isAcceptor[i] ) verlet->insert( i, &coords[i] );
  }
  verlet->makeNewList();
}

  HBManager::~HBManager() {
  delete firstHbond;
  delete verletCutoffs;
  delete verletCollector;
  delete verlet;
}

void HBManager::tighten() {

  for ( BBHBContainer::iterator bbhb = bb_hbonds.begin(); bbhb != bb_hbonds.end(); bbhb++ ) {
    tightenConstraint( *bbhb );
  }

  for ( HBContainer::iterator hb = hbonds.begin(); hb != hbonds.end(); hb++ ) {
    tightenConstraint( *hb );
  }
}

void HBManager::findnew( bool breakable ) {
  verlet->update();

  VerletCollectorHAPairs::const_iterator it;
  for ( it = verletCollector->begin(); it != verletCollector->end(); it++ ) {
    bool success;
    addConstraint( it->h, it->a, success, breakable );
  }
}

bool HBManager::isPairExcluded( int p1, int p2 ) const {
  int N = nt.size();
  if ( p1 >= N || p2 >= N ) return true;

  bool isPairExcluded;
  double cutoff;
  verletCutoffs->getCutoff( p1, p2, isPairExcluded, cutoff );
  return isPairExcluded;
}


void HBManager::addConstraint( int h, int a, bool& success, bool breakable ) {
  success = false;

  if ( forbidlist && forbidlist->checkforbid( h, a ) ) return;
  if ( nt[h].size() == 0 ) return;
  int d = nt[h][0];
  if ( forbidlist && forbidlist->checkforbid( d, a ) ) return;

  //check to see if this hbond already exists in our lists,
  //if so, don't add it again.
  bool isBB_Hbond = isBB[h] && isBB[a];
  if ( isBB_Hbond ) {
    if ( bb_hbonds.find( h, a ) != bb_hbonds.end() ) return;
  }
  else if ( hbonds.find( h, a ) != hbonds.end() ) return;

  if ( carefulChecking ) {
    //make sure energy is better than some minimal weak value.
    //This also ensures that the distance and angle meet
    //are within required parameters.
    double energy = firstHbond->energy( h, a );
    if ( energy > Ecutoff ) return;

    //Only one hydrogen bond is allowed per h atom.  If another h bond already
    //exists for this h atom, check to see if this new one is lower in energy.
    //if so, delete the old hbond and add the new hbond.
    vector<int> acceptors;
    bb_hbonds.lookupConstraintPartners( h, acceptors );
    if ( acceptors.size() > 0 ) {
      //if we find another hydrogen bond already exists, decide whether
      //to erase the old one and replace it with a new one,
      //or skip this new bond because it's not as good as the one that
      //already exists.
      BBHBContainer::iterator it = bb_hbonds.find( h, acceptors[0] );
      if ( it->isBreakable() && firstHbond->energy( h, acceptors[0]) > energy ) {
        bb_hbonds.erase( it );
        nRemoved_++;
      }
      else return;
    }
    //if we reach here, then any bb hydrogen bond that existed for this h
    //has been deleted.  But there still could be a regular hbond ( non-bb ).

    hbonds.lookupConstraintPartners( h, acceptors );
    if ( acceptors.size() > 0 ) {
      HBContainer::iterator it = hbonds.find( h, acceptors[0] );
      if ( it->isBreakable() && firstHbond->energy( h, acceptors[0]) > energy ) {
        hbonds.erase( it );
        nRemoved_++;
      }
      else return;
    }
    //if we reach here, then any hydrogen bond that may have existed for this h
    //has been deleted.
  }

  //We can now add the new hydrogen bond constraint.
  HBConstraint hb;
  BBHBConstraint bbhb;
  nAdded_++;
  if ( isBB_Hbond ) {
    initializeConstraint( h, a, bbhb );
    if ( !breakable ) bbhb.makeUnbreakable();
    bb_hbonds.insert( h, a, bbhb );
  }
  else {
    initializeConstraint( h, a, hb );
    if ( !breakable ) hb.makeUnbreakable();
    hbonds.insert( h, a, hb );
  }
  success = true;
}

void HBManager::initializeConstraint( int h, int a, HBConstraint& hb ) {
  //initialize atoms
  int d = nt[h][0];
  hb.setAtoms( &coords, d, h, a );

  //initialize dist
  hb.setMaxDistHA( max( L, hb.calcDistHA() ) );
  hb.setkDistHA( k );

  if ( doAngleConstraints )  {
    //initialize angle DHA
    hb.setMinAngleRadDHA( theta1 );
  }
  else {
    hb.setMinAngleRadDHA( 0 );
  }

  hb.setkAngleDHA( k );

}

void HBManager::initializeConstraint( int h, int a, BBHBConstraint& bbhb ) {
  //initialize atoms
  int d = nt[h][0];
  int b = nt[a][0];
  bbhb.setAtoms( &coords, d, h, a, b );

  //initialize dist
  bbhb.setMaxDistHA( max( L, bbhb.calcDistHA() ) );
  bbhb.setkDistHA( k );

  if ( doAngleConstraints )  {
    //initialize angle DHA
    double theta = bbhb.calcAngleRadDHA();
    double constraintAngle = theta >= theta2 ? theta2 : theta1;
    bbhb.setMinAngleRadDHA( constraintAngle );

    //initialize angle HAB
    theta = bbhb.calcAngleRadHAB();
    constraintAngle = theta >= 130 ? 130 : theta1;
    bbhb.setMinAngleRadHAB( constraintAngle );
  }
  else {
    bbhb.setMinAngleRadDHA( 0 );
    bbhb.setMinAngleRadHAB( 0 );
  }

  bbhb.setkAngleDHA( k );
  bbhb.setkAngleHAB( k );

}

void HBManager::tightenConstraint( BBHBConstraint& bbhb ) {
  if ( !bbhb.isBreakable() ) return;
  //tighten the distance constraint, but not tighter than L
  double constraintDist = bbhb.getConstraintMaxDistHA();
  if ( constraintDist > L + numeric_limits<double>::epsilon() ) {
    double dist = bbhb.calcDistHA();
    if ( dist < L ) bbhb.setMaxDistHA( L );
    else if ( dist < constraintDist - 0.1 ) {
      bbhb.setMaxDistHA( dist + 0.1 );
    }
  }

  if ( !doAngleConstraints ) return;

  //If the constraint is theta1, and if ever the DHA angle becomes > theta2,
  //update the constraint to theta2.
  double constraintAngle = bbhb.getConstraintMinAngleRadDHA();
  if ( constraintAngle < theta2 - numeric_limits<double>::epsilon() &&
       bbhb.calcAngleRadDHA() > theta2 ) {
    bbhb.setMinAngleRadDHA( theta2 );
  }

  //The HAB angle constraint is set to 0 at first.
  //As soon as the HAB angle becomes > theta2,
  //enable the constraint by setting it to theta2.
  if ( bbhb.getConstraintMinAngleRadDHA() < 130 - numeric_limits<double>::epsilon() &&
       bbhb.calcAngleRadHAB() >= 130 ) {
    bbhb.setMinAngleRadHAB( 130 );
  }
}

void HBManager::tightenConstraint( HBConstraint& hb ) {
  if ( !hb.isBreakable() ) return;

  //tighten the distance constraint, but not tighter than L
  double constraintDist = hb.getConstraintMaxDistHA();
  if ( constraintDist > L + numeric_limits<double>::epsilon() ) {
    double dist = hb.calcDistHA();
    if ( dist < L ) hb.setMaxDistHA( L );
    else if ( dist < constraintDist - 0.1 ) {
      hb.setMaxDistHA( dist + 0.1 );
    }
  }

  //we leave the angle constraint as is--no tightening
}

void HBManager::breakAllBreakable() {

  //first carry out for the BBHB constraints
  BBHBContainer::iterator it = bb_hbonds.begin();
  while ( it != bb_hbonds.end() ) {
    if ( it->isBreakable() ) {
      bb_hbonds.erase( it );
      nRemoved_++;
      //DO NOT advance the iterator "it".  The erase operation
      //of DynamicConstraintContainer replaces the current constraint
      //with the final constraint.  So, the "next" constraint
      //to be checked is the one that "it" is currently pointing to.
    }
    else it++;
  }

  //repeat for the HB constraints
  HBContainer::iterator ithb = hbonds.begin();
  while ( ithb != hbonds.end() ) {
    if ( ithb->isBreakable() ) {
      hbonds.erase( ithb );
      nRemoved_++;
    }
    else ithb++;
  }
}


/*
void DynamicHbonds::snipWeakBonds( double e ) {
  for ( int i = 0; i < nropes; i++ ) {
    if ( !(*ropes)[i].isEnabled ) continue;
    int p1 = (*ropes)[i].p1;
    int p2 = (*ropes)[i].p2;
    if ( isDonorH[p1] ) {
      double energy = firsthbond->energy( p1, p2 );
      if ( energy > -1.0 ) {
  (*ropes)[i].snip();
        map< pair<int,int>, double >::iterator it;
        it = hbondLookup.find( pair<int,int>( p1, p2 ) );
        hbondLookup.erase( it );
  continue;
      }
    }
    else if ( isDonorH[p2] ) {
      double energy = firsthbond->energy( p2, p1 );
      if ( energy > -1.0 ) {
  (*ropes)[i].snip();
        map< pair<int,int>, double >::iterator it;
        it = hbondLookup.find( pair<int,int>( p1, p2 ) );
        hbondLookup.erase( it );
  continue;
      }
    }
  }
}

void DynamicHbonds::snipMask( const vector<char>& mask ) {
  for ( int i = 0; i < nropes; i++ ) {

    if ( solvExp[p1] && solvExp[p2] ) {
      (*ropes)[i].snip();
      if ( isDonorH[p1] ) {
        map< pair<int,int>, double >::iterator it;
        it = hbondLookup.find( pair<int,int>( p1, p2 ) );
        hbondLookup.erase( it );
      }
      else if ( isDonorH[p2] ) {
        map< pair<int,int>, double >::iterator it;
        it = hbondLookup.find( pair<int,int>( p2, p1 ) );
        hbondLookup.erase( it );
      }
      else {
        set< pair<int,int> >::iterator it;
        it = phobeLookup.find( pair<int,int>( p1, p2 ) );
        phobeLookup.erase( it );
      }
      outfile << p1 << " " << p2 << " 0 " << confNumber << '\n';
      countSnipped++;
    }
  }
}
*/
