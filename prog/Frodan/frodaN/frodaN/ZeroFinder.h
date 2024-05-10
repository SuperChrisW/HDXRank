#ifndef ZEROFINDER_H_
#define ZEROFINDER_H_

#include <cmath>

template <class F, class Real >
class ZeroFinder {
 public:
  Real findZero( F &func, Real bracket1, Real bracket2, Real xtol );
};

template <class F, class Real >
Real ZeroFinder<F, Real >::findZero( F &func, Real x1, Real x2, Real xtol ) {
  //this function implements the bisection method of finding a zero.
  //It requires that the zero lies within the bracket [x1, x2],
  //and that the zero occurs by crossing through the x-axis
  //(not curving down to touch the x-axis without crossing)

  Real f1 = func(x1);
  Real f2 = func(x2);
  Real fmid;
  Real xmid;
  Real dx;

  
  //probably won't happen, but if one of the input bracket values
  //happens to be the zero, then we already have the answer
  if ( f1 == 0 ) return x1;
  else if ( f2 == 0 ) return x2;

  //ensure that bracket x1 gives the negative function value,
  //and bracket x2 gives the positive function value,
  //interchanging brackets if needed
  else if ( f1 > 0 && f2 < 0 ) {
    double temp = x2;
    x2 = x1;
    x1 = temp;
  }
  dx = x2-x1;

  //bisection method,
  //find the zero by successive bisections
  //of the bracket
  for ( int count = 0; count < 40; count++ ) {
    //calculate bracket midpoint
    dx *= 0.5;
    xmid = x1 + dx;

    //if the zero is sure to be found within + or - xtol of the midpoint, we are done
    if ( fabs(dx) < xtol ) break;

    //evaluate function at midpoint
    fmid = func(xmid);

    //if we found the zero exactly (probably a rare occurrence), we are done
    if ( fmid==0.0 ) break;

    //tighten bracket
    if ( fmid < 0 ) x1=xmid;
    else x2=xmid;
  }

  return xmid;
}



#endif
