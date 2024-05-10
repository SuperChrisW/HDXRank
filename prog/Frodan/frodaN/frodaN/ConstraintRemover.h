/*
 * ConstraintRemover.h
 *
 *  Created on: Dec 8, 2009
 *      Author: dwfarrel
 */

#ifndef CONSTRAINTREMOVER_H_
#define CONSTRAINTREMOVER_H_

#include "ConstraintEnforcingPotential_Targeting.h"
#include <queue>
#include "ForbidList.h"

class ConstraintRemoverAbstract {
public:
  ConstraintRemoverAbstract( ConstraintEnforcingPotential_Targeting* cep_ ) :
    cep( cep_ ), forbidlist( NULL ) {}
  virtual ~ConstraintRemoverAbstract() {}
  virtual void remove( int& nRemoved ) = 0;
  void attachForbidList( ForbidList* forbidlist_ ) {
    forbidlist = forbidlist_;
  }
protected:
  ConstraintEnforcingPotential_Targeting* cep;
  ForbidList* forbidlist;
};

class ConstraintRemover_Random : public ConstraintRemoverAbstract {
public:
  ConstraintRemover_Random( 
    ConstraintEnforcingPotential_Targeting* cep_ ) :
      ConstraintRemoverAbstract( cep_ ),
      f_bbhb( 0 ),
      f_hb( 0 ),
      f_ph( 0 ) {}
    void setRemoveFrac_bbhb( double f ) { f_bbhb = f; }
    void setRemoveFrac_hb( double f ) { f_hb = f; }
    void setRemoveFrac_ph( double f ) { f_ph = f; }
  virtual ~ConstraintRemover_Random() {}
  void remove( int& nRemoved );
private:
  double f_bbhb;
  double f_hb;
  double f_ph;
};

class ConstraintRemover_Cutoff : public ConstraintRemoverAbstract {
public:
  ConstraintRemover_Cutoff( ConstraintEnforcingPotential_Targeting* cep_ ) :
    ConstraintRemoverAbstract( cep_ ), d_bbhb( 0 ), d_hb( 0 ), d_ph( 0 ) {}
  virtual ~ConstraintRemover_Cutoff() {}
  void setStretchCutoff_bbhb( double d_bbhb_ ) { d_bbhb = d_bbhb_; }
  void setStretchCutoff_hb( double d_hb_ ) { d_hb = d_hb_; }
  void setStretchCutoff_ph( double d_ph_ ) { d_ph = d_ph_; }
  void remove( int& nRemoved );
private:
  double d_bbhb;
  double d_hb;
  double d_ph;
};

#endif /* CONSTRAINTREMOVER_H_ */
