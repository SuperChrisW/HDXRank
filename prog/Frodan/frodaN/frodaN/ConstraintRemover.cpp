/*
 * ConstraintRemover.cpp
 *
 *  Created on: Dec 8, 2009
 *      Author: dwfarrel
 */

#include "ConstraintRemover.h"
#include "mt19937ar.h"

using namespace std;

void ConstraintRemover_Random::remove( int& nRemoved ) {
  //check each constraint.
  //if it is breakable,
  //remove it with the pre-set probability.
  nRemoved = 0;

  if ( cep->bbhb ) {
    BBHBContainer* bbhb = cep->bbhb;
    BBHBContainer::iterator itbbhb = bbhb->begin();
    while ( itbbhb != bbhb->end() ) {
      int p1 = itbbhb->geth();
      int p2 = itbbhb->geta();
      int d = itbbhb->getd();
      bool isBreakable = itbbhb->isBreakable();

      if ( isBreakable && genrand_real2() < f_bbhb ) {
        if ( forbidlist ) forbidlist->inserthb( d, p2 );
        bbhb->erase( itbbhb );
        //DO NOT increment itbbhb if we are erasing, because the constraint container
        //automatically moves a constraint into the vacant position.
        nRemoved++;
      }
      else itbbhb++;
    }
  }

  if ( cep->hb ) {
    HBContainer* hb = cep->hb;
    HBContainer::iterator ithb = hb->begin();
    while ( ithb != hb->end() ) {
      int p1 = ithb->geth();
      int p2 = ithb->geta();
      int d = ithb->getd();
      bool isBreakable = ithb->isBreakable();

      if ( isBreakable && genrand_real2() < f_hb ) {
        if ( forbidlist ) { forbidlist->inserthb( d, p2 ); }
        hb->erase( ithb );
        nRemoved++;
      }
      else ithb++;
    }
  }

  if ( cep->ph ) {
    PHContainer* ph = cep->ph;
    PHContainer::iterator itph = ph->begin();
    while ( itph != ph->end() ) {
      int p1 = itph->getp1();
      int p2 = itph->getp2();
      bool isBreakable = itph->isBreakable();

      if ( isBreakable && genrand_real2() < f_ph ) {
        if ( forbidlist ) forbidlist->insert( p1, p2 );
        ph->erase( itph );
        //DO NOT increment itbbhb if we are erasing, because the constraint container
        //automatically moves a constraint into the vacant position.
        nRemoved++;
      }
      else itph++;
    }
  }

  if ( nRemoved ) { cep->notifyTermsChanged(); }
}


void ConstraintRemover_Cutoff::remove( int& nRemoved ) {
  //check each constraint.
  //if it is breakable,
  //and if it is over-stretched by some threshold level,
  //remove it with the pre-set probability.
  nRemoved = 0;
  if ( cep->bbhb ) {
    BBHBContainer* bbhb = cep->bbhb;
    BBHBContainer::iterator itbbhb = bbhb->begin();
    while ( itbbhb != bbhb->end() ) {
      if ( itbbhb->isBreakable() &&
           itbbhb->calcDistHA() - itbbhb->getConstraintMaxDistHA() > d_bbhb ) {
        if ( forbidlist ) forbidlist->inserthb( itbbhb->getd(), itbbhb->geta() );
        bbhb->erase( itbbhb );
        nRemoved++;
      }
      else itbbhb++;
    }
  }

  if ( cep->hb ) {
    HBContainer* hb = cep->hb;
    HBContainer::iterator ithb = hb->begin();
    while ( ithb != hb->end() ) {
      if ( ithb->isBreakable() &&
           ithb->calcDistHA() - ithb->getConstraintMaxDistHA() > d_hb ) {
        if ( forbidlist ) forbidlist->inserthb( ithb->getd(), ithb->geta() );
        hb->erase( ithb );
        nRemoved++;
      }
      else ithb++;
    }
  }

  if ( cep->ph ) {
    PHContainer* ph = cep->ph;
    PHContainer::iterator itph = ph->begin();
    while ( itph != ph->end() ) {
      if ( itph->isBreakable() &&
           itph->calcDist() - itph->getCutoff() > d_ph ) {
        if ( forbidlist ) forbidlist->insert( itph->getp1(), itph->getp2() );
        ph->erase( itph );
        nRemoved++;
      }
      else itph++;
    }
  }

  if ( nRemoved ) { cep->notifyTermsChanged(); }
}
