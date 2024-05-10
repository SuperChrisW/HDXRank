/*
 * CollectiveRandomPerturber.cpp
 *
 *  Created on: Apr 27, 2010
 *      Author: dwfarrel
 */

#include "CollectiveRandomPerturber.h"
#include "mt19937ar.h"
#include "RandomVector.h"
#include "ProteinInfo.h"
#include "NeighborTable.h"
#include "RigidUnitSystem.h"
#include "ConstraintEnforcingPotential.h"
#include <set>
#include <queue>
#include <cmath>

using namespace std;

CollectiveRandomPerturber::CollectiveRandomPerturber(
  const ProteinInfo* prot_,
  const NeighborTable* covNT_,
  RigidUnitSystem* sys_ ,
  const ConstraintEnforcingPotential* cep_ ) :
  prot( prot_ ),
  covNT( *covNT_ ),
  sys( sys_ ),
  cep( cep_ ),
  resNT( prot->nresi() )
{

  //build a "neighbor table" for residues, that connects each residue
  //with other residues that are covalently bonded
  //OR backbone hydrogen bonded (add that part later)
  int natoms = prot->natoms();
  for ( int atom1 = 0; atom1 < natoms; atom1++ ) {
    int res1 = prot->atom(atom1).resi().index();
    int Nneigh = covNT[atom1].size();
    for ( int j = 0; j < Nneigh; j++ ) {
      int atom2 = covNT[atom1][j];
      if ( atom1 >= atom2 ) continue; //only check pairs once (i<j)
      int res2 = prot->atom(atom2).resi().index();
      if ( res1 != res2 ) {
        resNT.insert( res1, res2 );
      }
    }
  }

  if ( cep->bbhb ) {
    BBHBContainer* bbhb = cep->bbhb;
    for ( BBHBContainer::iterator itbbhb = bbhb->begin(); itbbhb != bbhb->end(); itbbhb++ ) {
      if ( itbbhb->isBreakable() ) continue;
      int res1 = prot->atom( itbbhb->geth() ).resi().index();
      int res2 = prot->atom( itbbhb->geta() ).resi().index();
      if ( res1 != res2 ) {
        resNT.insert( res1, res2 );
      }
    }
  }

  if ( cep->hb ) {
    HBContainer* hb = cep->hb;
    for ( HBContainer::iterator ithb = hb->begin(); ithb != hb->end(); ithb++ ) {
      if ( ithb->isBreakable() ) continue;
      int res1 = prot->atom( ithb->geth() ).resi().index();
      int res2 = prot->atom( ithb->geta() ).resi().index();
      if ( res1 != res2 ) {
        resNT.insert( res1, res2 );
      }
    }
  }

  if ( cep->ph ) {
    PHContainer* ph = cep->ph;
    for ( PHContainer::iterator itph = ph->begin(); itph != ph->end(); itph++ ) {
      if ( itph->isBreakable() ) continue;
      int res1 = prot->atom( itph->getp1() ).resi().index();
      int res2 = prot->atom( itph->getp2() ).resi().index();
      if ( res1 != res2 ) {
        resNT.insert( res1, res2 );
      }
    }
  }


  resNT.commit();

  //build RUlistFromResi
  int nresi = prot->nresi();
  RUlistFromResi.resize( nresi );
  for ( int r = 0; r < nresi; r++ ) {
    set<int> ru_set;
    ru_set.clear();
    Resi::const_iterator resbegin = prot->resi(r).begin();
    Resi::const_iterator resend = prot->resi(r).end();
    for ( Resi::const_iterator atom = resbegin; atom != resend; atom++ ) {
      //for each residue, we only want to associate the rigid units of the side chain
      //and C-alpha and the peptide plane of the N-side.  We omit the peptide plane
      //of the C-side.
      if ( atom->name() == "C" || atom->name() == "O" ) continue;

      int atomindex = atom->index();
      ru_set.insert( sys->getRUlistFromP( atomindex ).begin(), sys->getRUlistFromP( atomindex ).end() );
    }
    RUlistFromResi[r].insert( RUlistFromResi[r].begin(), ru_set.begin(), ru_set.end() );
  }

}

CollectiveRandomPerturber::~CollectiveRandomPerturber() {
}

void CollectiveRandomPerturber::perturb() {
  //initialize bucket of residues,
  //and initialize random ordering for picking residues.
  resi_bucket.clear();
  priority_queue<PQNode> pq;
  int Nresi = prot->nresi();
  for ( int i = 0; i < Nresi; i++ ) {
    resi_bucket.insert( i );
    pq.push( PQNode( genrand_real2(), i ) );
  }

  //pick out one residue randomly and all neighbor residues out to xth neighbor
  while ( !pq.empty() ) {
    vector<int> resi_group;

    int resi = pq.top().resi;
    pq.pop();

    resi_group.clear();
    extendGroup( resi, 0, resi_group );

    //convert resi group to a rigid unit group
    set<int> rigidUnitGroup;
    int N = resi_group.size();
    for ( int i = 0; i < N; i++ ) {
      int resi = resi_group[i];
      rigidUnitGroup.insert( RUlistFromResi[resi].begin(), RUlistFromResi[resi].end() );
    }

    //perturb the chosen group of residues
    if ( rigidUnitGroup.size() == 0 ) continue;
    perturbRigidUnitGroup_Rotational( rigidUnitGroup );
    //possibly the rotational pert has to come first,
    //since it requires that the mean atom positions be updated/in-sync with the
    //embedded rigid unit positions.
    //The translational perturbation does not have this requirement.
    perturbRigidUnitGroup_Translational( rigidUnitGroup );

  }

  sys->update();

}

void CollectiveRandomPerturber::extendGroup( int resi, int nExtend, vector<int>& resi_group ) {
  //recursive function

  //if resi is not in the bucket, stop
  set<int>::iterator it = resi_bucket.find( resi );
  if ( it == resi_bucket.end() ) return;

  //transfer from bucket to resi_group
  resi_group.push_back( resi );
  resi_bucket.erase( it );

  //if no more extensions, stop
  if ( nExtend <= 0 ) return;

  //extend group along each neighbor
  int nNeigh = resNT[resi].size();
  for ( int k = 0; k < nNeigh; k++ ) {
    int neighres = resNT[resi][k];
    extendGroup( neighres, nExtend-1, resi_group );
  }

}

void CollectiveRandomPerturber::perturbRigidUnitGroup_Translational( const set<int>& rigidUnitGroup ) {
  Vec3 centerPerturbation;
  double pertTranslationSize = 1.0;
  generateRandomUnitVector( centerPerturbation );
  centerPerturbation *= genrand_real2()*pertTranslationSize;

  for ( set<int>::iterator it = rigidUnitGroup.begin(); it != rigidUnitGroup.end(); it++ ) {
    sys->addToCenter( *it, centerPerturbation );
  }

}

void CollectiveRandomPerturber::perturbRigidUnitGroup_Rotational( const set<int>& rigidUnitGroup ) {
  double arcLengthPertSize = 1.0;

  Vec3 rotAxis;
  generateRandomUnitVector( rotAxis );

  set<int> atom_set;

  for ( set<int>::iterator it = rigidUnitGroup.begin(); it != rigidUnitGroup.end(); it++ ) {
    int ru = *it;
    atom_set.insert( sys->getPlistFromRU( ru ).begin(), sys->getPlistFromRU( ru ).end() );
  }

  Vec3 center(0,0,0);
  for ( set<int>::iterator it = atom_set.begin(); it != atom_set.end(); it++ ) {
    center += sys->meanPositions( *it );
  }
  center /= atom_set.size();

  double maxProjectionPerp2 = 0;
  for ( set<int>::iterator it = atom_set.begin(); it != atom_set.end(); it++ ) {
    Vec3 rel = sys->meanPositions( *it ) - center;
    double mag2 = rel.norm2();
    double projection_alongAxis = rel.dot( rotAxis );
    double projection_perp2 = mag2 - projection_alongAxis*projection_alongAxis;
    if ( projection_perp2 > maxProjectionPerp2 ) maxProjectionPerp2 = projection_perp2;
  }
  double maxProjectionPerp = sqrt( maxProjectionPerp2 );


  double maxRotationAngle = arcLengthPertSize/maxProjectionPerp;
  double randomRotationAngle = genrand_real1()*maxRotationAngle;
  double rotorMagnitude = 2.0 * sin(randomRotationAngle/2.0);

  Rotator rotator( rotAxis*rotorMagnitude );
  for ( set<int>::iterator it = rigidUnitGroup.begin(); it != rigidUnitGroup.end(); it++ ) {
    int ru = *it;
    sys->rotate( ru, rotator, center );
  }

}
