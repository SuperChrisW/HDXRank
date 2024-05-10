/*
 * OverlapEnergy.cpp
 *
 *  Created on: Jul 24, 2009
 *      Author: dwfarrel
 */

#include "Rama.h"
#include "ProteinInfo.h"
#include "NeighborTable.h"
#include "RigidUnitSystem.h"
#include <cmath>
#include <limits>

using namespace std;

RamaContainer::RamaContainer(
    const RigidUnitSystem* sys_,
    const ProteinInfo& proteinInfo_,
    const NeighborTable& nt_ ) :
  ConstraintContainer<MinDistConstraint>( proteinInfo_.natoms() ),
  sys(sys_),
  proteinInfo(proteinInfo_),
  nt( nt_ ),
  k( 10.0 )
  //k( 20.0 )
{
  //only consider atoms within the standard amino acids,
  //GLY is included in this list, but because it lacks a CB atom,
  //only the Rama constraints for the atoms it does have will
  //will be established.
  resList.insert( "LEU" );
  resList.insert( "ALA" );
  resList.insert( "VAL" );
  resList.insert( "MET" );
  resList.insert( "ILE" );
  resList.insert( "TRP" );
  resList.insert( "PHE" );
  resList.insert( "SER" );
  resList.insert( "THR" );
  resList.insert( "CYS" );
  resList.insert( "ASN" );
  resList.insert( "GLN" );
  resList.insert( "TYR" );
  resList.insert( "HIS" );
  resList.insert( "ASP" );
  resList.insert( "GLU" );
  resList.insert( "LYS" );
  resList.insert( "ARG" );
  resList.insert( "GLY" );
  resList.insert( "PRO" );

  //constraints producing rectangular shaped phi disallowed region
  double constraint_Oprev_CB = 3.05; //suggested increase from 3.05 to 3.1

  //constraints producing phi-psi upper oval disallowed region
  double constraint_Oprev_O = 3.1; //suggested increase from 3.1 to 3.3
  double constraint_Cprev_O = 2.65; //taken as 95% of pair type cutoff for this pair
  //constraints producing phi-psi central oval disallowed region
  double constraint_Oprev2_N = 3.00; //suggested increase from 3.0 to 3.2
  //constraint producing phi-psi left-side small oval
  double constraint_Hprev_H = 1.85; //suggested increase from 1.85 to 1.9
  //double constraint_Hprev_N = 2.50; //NEW 2.5

  //constraint producing psi disfavored rectangular region 1
  double constraint_CB_O = 2.8; //increase from 2.8 to 3.0
  //constraint producing psi disfavored rectangular region 2
  double constraint_CBprev_N = 3.0; //increase from 3.0 to 3.1


  //for each residue, lookup indices for mainchain atoms of interest
  //verify a peptide bond between residues
  int cb = -1;
  int c = -1;
  int o = -1;
  int n = -1;
  int h = -1;

  int cbprev = -1;
  int cprev = -1;
  int oprev = -1;
  int nprev = -1;
  int hprev = -1;

  int oprev2 = -1;

  bool peptidebond_prev_to_prev2 = false;
  bool peptidebond_current_to_prev = false;

  for ( ProteinInfo::Residues::const_iterator r = proteinInfo.residues().begin();
        r != proteinInfo.residues().end(); r++ ) {
    oprev2 = oprev;

    cbprev = cb;
    cprev = c;
    oprev = o;
    nprev = n;
    hprev = h;

    cb = -1;
    c = -1;
    o = -1;
    n = -1;
    h = -1;
    for ( Resi::const_iterator atom = r->begin(); atom != r->end(); atom++ ) {
      if ( atom->name() == "CB" ) cb = atom->index();
      else if ( atom->name() == "C" ) c = atom->index();
      else if ( atom->name() == "O" ) o = atom->index();
      else if ( atom->name() == "N" ) n = atom->index();
      else if ( atom->name() == "H" || atom->name() == "HN" ) 
        h = atom->index();
    }

    peptidebond_prev_to_prev2 = peptidebond_current_to_prev;
    peptidebond_current_to_prev = cprev >= 0 && n >= 0 && nt.isFirstNeighbor( cprev, n );

    if ( cb >= 0 && o >= 0 ) {
      attemptAddRamaConstraint( cb, o, constraint_CB_O );
    }
    if ( peptidebond_current_to_prev ) {
      if ( cbprev >= 0 && n >= 0 ) attemptAddRamaConstraint( cbprev, n, constraint_CBprev_N );
      if ( oprev >= 0 && o >= 0 ) attemptAddRamaConstraint( oprev, o, constraint_Oprev_O );
      if ( cprev >= 0 && o >= 0 ) attemptAddRamaConstraint( cprev, o, constraint_Cprev_O );
      if ( oprev >= 0 && cb >= 0 ) attemptAddRamaConstraint( oprev, cb, constraint_Oprev_CB );
      if ( hprev >= 0 && h >= 0 ) attemptAddRamaConstraint( hprev, h, constraint_Hprev_H );
      //if ( hprev >= 0 && n >= 0 ) attemptAddRamaConstraint( hprev, n, constraint_Hprev_N );
      if ( peptidebond_prev_to_prev2 ) {
        if ( oprev2 >= 0 && n >= 0 ) attemptAddRamaConstraint( oprev2, n, constraint_Oprev2_N );
      }
    }
  }
}

RamaContainer::~RamaContainer() {
}

void RamaContainer::attemptAddRamaConstraint( int i, int j, double c ) {
  if ( i == j ||
       sys->doPointsBelongToSameRigidUnit( i, j ) ||
       nt.isFirstNeighbor( i, j ) ||
       nt.isSecondNeighbor( i, j ) ||
       resList.find( proteinInfo.atom(i).resi().name() ) == resList.end() ||
       resList.find( proteinInfo.atom(j).resi().name() ) == resList.end() ) return;
  insert( i, j, MinDistConstraint( &sys->meanPositions(), k, i, j, c ) );
}
