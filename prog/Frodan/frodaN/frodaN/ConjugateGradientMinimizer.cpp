#include "ConjugateGradientMinimizer.h"
#include "MultiVarFunction.h"
#include "MultiVarFunction_1DProjection.h"
#include <limits>
#include <iomanip>
#include <algorithm>
#include <iostream>

using namespace std;

ConjugateGradientMinimizer::ConjugateGradientMinimizer() :
  f(NULL),
  lambda( 0 ),
  reachedTol( false ),
  step(0),
  isEnabled_TolMaxGradComp( false ),
  isEnabled_TolMaxPreconditionedGradComp( false ),
  tolMaxGradComp( 0.0 ),
  tolMaxPreconditionedGradComp( 0.0 ),
  outputPeriod( 0 ),
  doSteepestDescent_or_doConjGrad_0or1( 1 ),
  requireRigorousLineMinimization( true ),
  trialQstep( 0.01 ),
  maxQstep( 5 ),
  flag_foundZeroGradient( false )
{
}

ConjugateGradientMinimizer::~ConjugateGradientMinimizer()
{
}

bool ConjugateGradientMinimizer::toleranceSatisfied() {
  if ( flag_foundZeroGradient ) return true;
  if ( isEnabled_TolMaxGradComp ) {
    double maxgradcomp = 0.0;
    size_t N = r.size();
    for ( size_t i=0; i<N; i++ ) {
      double ri = abs(r[i]);
      if ( ri > maxgradcomp ) { maxgradcomp = ri; }
    }
    if ( maxgradcomp > tolMaxGradComp ) return false;
  }
  if ( isEnabled_TolMaxPreconditionedGradComp ) {
    if ( !isPreconditioned ) return false;
    double maxs = 0.0;
    size_t N = s.size();
    for ( size_t i=0; i<N; i++ ) {
      double si = abs(s[i]);
      if ( si > maxs ) { maxs = si; }
    }
    if ( maxs > tolMaxPreconditionedGradComp ) return false;
  }
  if ( !isEnabled_TolMaxGradComp && !isEnabled_TolMaxPreconditionedGradComp ) {
    cout << "Error: No tolerance condition set for minimizer" << endl;
    exit(0);
  }
  return true;
}

void ConjugateGradientMinimizer::stepAlongCurrentDirection()
{
  //assume that the next direction p is already set

  //cout << "stepAlongCurrentDirection Begin: " << "\n";

  //find pmax, the largest component in the new direction p.
  size_t N = p.size();
  double pmax = 0.0;
  for ( size_t i=0; i<N; i++ ) {
    double absp = abs(p[i]);
    if ( absp > pmax ) { pmax = absp; }
  }

  if ( pmax*pmax < numeric_limits<double>::epsilon() ) {
    lambda = 0;
    flag_foundZeroGradient = true;
    return;
  }

  //Note that lambda is the step size, multiplying direction p.
  //The step vector then is q = lambda*p.
  //If the largest component of vector p is some number P, then
  //The largest component Q of the step vector q is Q = lambda*P.
  //Hence, lambda/Q = 1/P.
  double lambdaOverQ = 1.0/pmax ;

  //to set the first trial step size lambda,
  //consider how large a step in the largest component Q we can tolerate on the trial step.
  //This number is trialQstep.
  //The lambda that will give us this step in component Q is given by
  //trialQstep*LambdaOverQ.
  double trialLambda = trialQstep*lambdaOverQ;
  double lambdaMax = maxQstep*lambdaOverQ; //1000.0;

  //double trialLambda = min(trial_dqmax*dLambda_dQ, trial_lambdamax );
  findBracket( multiVarFunction_1DProjection, 0.0, trialLambda, lambdaMax );
  //cout << "FindBracket end" << endl;

  if ( findBracket.allowEarlyStop &&
        findBracket.likelyConvergedToMinimum ) {
    if ( findBracket.foundBracket ) lambda = findBracket.b;
    else lambda = findBracket.c;
  }
  else if ( findBracket.foundBracket ) {
    if ( findBracket.a > findBracket.c ) {
      swap( findBracket.a, findBracket.c );
      swap( findBracket.fa, findBracket.fc );
    }

    //double tol = 1.0e-3 * findBracket.b;
    lambda = linemin(
                multiVarFunction_1DProjection, 100, 1.0e-2,
                findBracket.a, findBracket.b, findBracket.c,
                findBracket.fa, findBracket.fb, findBracket.fc );
    if ( !linemin.foundMinimum ) {
      cout << "Error: linemin failed " << endl;
      exit(0);
    }

    //cout << "LineMinimizer " << linemin.iteration << endl;
    //cout << "Minimum found.  Lambda = " << lambda << endl;
  }
  else if ( findBracket.foundFlatRegion ) {
    lambda = findBracket.b;
    cout << "Flat region found.  Lambda = " << lambda << endl;
  }
  else if ( findBracket.foundDecrease ) {
    lambda = lambdaMax;
    /*
    cout << "ConjGrad Error: Function only decreases along direction p," << endl;
    cout << "   Could not find a minimum, or a flat region" << endl;
    exit(0);
    */
  }
  else {
    lambda = 0.0;
    cout << "Error: Could not find a decrease along direction p" << endl;
    exit(0);
  }

  if ( lambda < 0.0 ) {
    /*
    cout << "Error: lambda is negative" << endl;
    exit(0);
    */
  }
  multiVarFunction_1DProjection.update( lambda );
}

void ConjugateGradientMinimizer::setNextDirection_SteepestDescent() {

  // set r to the new negative gradient
  f->getNegGrad( r );
  f->getPreconditionedNegGrad( s, isPreconditioned );

  p = s;
  multiVarFunction_1DProjection.setProjection( f, p );
}

void ConjugateGradientMinimizer::setNextDirection_ConjGrad() {
  // set r to the new negative gradient
  f->getNegGrad( r );
  r_dot_s_old = r_dot_s_new;

  r_dot_s_mid = 0;
  size_t N = r.size();
  for ( size_t i=0; i<N; i++ ) {
    r_dot_s_mid += r[i]*s[i];
  }

  f->getPreconditionedNegGrad( s, isPreconditioned );

  r_dot_s_new = 0.0;
  for ( size_t i=0; i<N; i++ ) {
    r_dot_s_new += r[i]*s[i];
  }

  double beta;
  if ( r_dot_s_old < numeric_limits<double>::epsilon() ) {
    beta = 0;
    flag_foundZeroGradient = true;
  }
  else {
    beta = (r_dot_s_new - r_dot_s_mid)/r_dot_s_old;
  }
  //cout << "ConjGrad beta " << beta << "\n";

  if ( beta > 0 && step_CGreset < 1000 ) {
    // p gets  (r + beta*p)
    for ( size_t i=0; i<N; i++ ) {
      p[i] = s[i] + beta*p[i];
    }
    step_CGreset++;
  }
  else {
    p = s;
    step_CGreset = 0;
  }

  // Now that the new p is calculated, set up the 1D model for the
  // next step
  multiVarFunction_1DProjection.setProjection( f, p );
}

void ConjugateGradientMinimizer::startMinimization( MultiVarFunction *f_ ) {
  f = f_;
  findBracket.allowEarlyStop = !requireRigorousLineMinimization;
  flag_foundZeroGradient = false;

  //to begin, r2 gets the current gradient's squared magnitude,
  //and direction p is set to be the negative gradient
  f->getNegGrad( r );
  f->getPreconditionedNegGrad( s, isPreconditioned );

  double maxr = 0;

  r_dot_s_new = 0.0;
  size_t N = r.size();
  for ( size_t i=0; i<N; i++ ) {
    double absr = abs( r[i] );
    if ( absr > maxr ) maxr = absr;
    r_dot_s_new += r[i]*s[i];
  }
  r_dot_s_old = r_dot_s_new;
  if ( maxr*maxr < numeric_limits<double>::epsilon() ) flag_foundZeroGradient = true;

  p = s;
  multiVarFunction_1DProjection.setProjection( f, p );

  reachedTol = toleranceSatisfied();

  step_CGreset=0;
  step=0;
}

void ConjugateGradientMinimizer::doStep() {
  stepAlongCurrentDirection();
  if ( doSteepestDescent_or_doConjGrad_0or1 == 0 ) {
    setNextDirection_SteepestDescent();
  }
  else {
    setNextDirection_ConjGrad();
  }
  reachedTol = toleranceSatisfied();
  step++;
}

void ConjugateGradientMinimizer::minimize( MultiVarFunction *f_, int nsteps )
{
  startMinimization( f_ );

  if ( outputPeriod ) {
    cout << "\nConjugate Gradient Minimization\n";
    cout << "  Nsteps = " << nsteps << '\n';
    cout << "  Reporting Every " << outputPeriod << " steps" << endl;
    cout << endl;
    cout << "Step | Lambda | Energy" << endl;
  }


  while ( !reachedTol && step<nsteps ) {
    doStep();
    if ( outputPeriod && ( step%outputPeriod == 0 ) ) {
      cout << step << " " << lambda << " " << f->eval() << endl;
    }
    notifyObservers();
  }

  if ( outputPeriod ) {
    if ( reachedTol ) {
      cout << "Tolerance Achieved after " << step << " steps." << endl;
    }
    else {
      cout << "Max number of steps exceeded (" << nsteps << ").  Tolerance Not Achieved." << endl;
    }
  }
}

