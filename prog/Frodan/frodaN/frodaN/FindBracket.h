#ifndef FINDBRACKET_H_
#define FINDBRACKET_H_

#include <iostream>
#include <algorithm>
#include <cmath>
#include <limits>

template <typename T>
class FindBracket
{
public:
	FindBracket() :
    allowEarlyStop( false ),
    beta(2 - 0.5*(3.0 - sqrt(5.0))), //beta gets 1.618
    limitingParabolicMagnification(10.0) {}
	virtual ~FindBracket() {}

  void operator()( T &f,
									 double a_,
									 double b_,
									 double xmax);

  /*
  void operator()( T &f,   double xinitial,
                    double deltax,
                    double beta,
                    double xmax,
                    double deltaf_flat );
                    */
  bool allowEarlyStop;
  double a;
  double b;
  double c;
  double fa;
  double fb;
  double fc;
  bool foundBracket;
  bool foundFlatRegion;
  bool foundDecrease;
  bool likelyConvergedToMinimum;
  int count;
  double beta;
  double limitingParabolicMagnification;
};

template <typename T>
void FindBracket<T>::operator()( T &f,
                                  double a_,
                                  double b_,
                                  double xmax) {
  count = 0;
  foundBracket = false;
  foundFlatRegion = false;
  foundDecrease = false;
  likelyConvergedToMinimum = false;
  a = a_;
  //cout << "Lambda a " << a << endl;
  fa = f(a);
  b = b_;
  //cout << "Lambda b " << b << endl;
  fb = f(b);
  if ( fb > fa ) {
    //return; //did not find an initial decrease
    std::swap(a,b);
    std::swap(fa,fb);
  }
  else {
    foundDecrease = true;
  }

  c = b + beta*(b-a);
  //cout << "Lambda c " << c << endl;
  fc = f(c);
  //we know the function starts out decreasing.
  //Now try to find an increase.
  while ( c < xmax ) {

    count++;
    double limit = b + limitingParabolicMagnification*(c-b);
    double fba = fb-fa;
    double fbc = fb-fc;

    if ( std::abs(fba) < std::numeric_limits<double>::epsilon() &&
         std::abs(fbc) < std::numeric_limits<double>::epsilon() ) {
      foundFlatRegion = true;
      return;
    }
    if ( fbc < 0.0 ) {
      foundBracket = true;
      return;
    }

    double ba = b-a;
    double bc = b-c;
    double u = b - 0.5*( bc*bc*fba - ba*ba*fbc )/( bc*fba - ba*fbc );
    double fu;
    //cout << " a " << a << " b " << b << " c " << c <<
      //   " u " << u <<
      //   " fa " << fa << " fb " << fb << " fc " << fc << endl;

    const double frac_tol = 0.001;
    if ( (( b<c && u>b ) || ( b>c && u<b )) &&
        std::abs((u-c)/c) < frac_tol ) {

      //u is very close to c (could be slightly greater or less than c).
      //This is fortunate, because although we are simply trying to bracket
      //the minimum, it is quite likely that we have actually
      //found the precise minimum.  If the boolean setting "allowEarlyStop"
      //is set to true, then we will immediately stop trying to bracket
      //the minimum.  The triple (a,b,c) does not bracket the minimum, but
      //the minimum is likely to be located at (or near) point c.  The
      //benefit of stopping immediately is that we save time, avoiding more
      //energy calculations to rigorously bracket the minimum when we
      //already are at the minimum.
      //(If that is too uncertain for you, then keep the allowEarlyStop
      //set to false!)
      if ( allowEarlyStop ) {
        likelyConvergedToMinimum = true;
        return;
      }

      //With u and c so close, we may not be able to distinguish their function
      //values fu and fc to tell if we have bracketed the minimum.
      //So, move c to be a little past u, enough that we can
      //distinguish the function values fc and fu, which will hopefully
      //allow us to succesfully bracket the minimum.
      c = (c>b) ? (1.0+frac_tol)*u : (1.0-frac_tol)*u;
      //cout << "Lambda c (just past u) " << c << endl;
      fc=f(c);
      //cout << " a " << a << " b " << b << " c " << c <<
        // " u " << u <<
        // " fa " << fa << " fb " << fb << " fc " << fc << endl;
    }

    if ( (b<c && u>b && u<c) ||
         (b>c && u<b && u>c) ) {
      //cout << "Lambda u (between b and c) " << u << endl;
      fu = f(u);
      if ( fu <= fc ) {
        //bracket found: (a,b,c) <- (b,u,c)
        a = b; fa = fb;
        b = u; fb = fu;
        foundBracket = true;
        //cout << " a " << a << " b " << b << " c " << c <<
          //" fa " << fa << " fb " << fb << " fc " << fc << endl;

        //before we leave the function, it is worth checking here
        //if our bracket b value is close to the parabolic reverse-interpolation
        //fit of points a,b,c.  If so, then we set the "likelyConvergedToMinimum"
        //flag to true, to inform the user.  (If we don't inform the user,
        //then he/she must assume that the minimum could be anywhere in
        //the bracket, between a and c).
        bc = b-c;
        ba = b-a;
        fbc = fb-fc;
        fba = fb-fa;
        u = b - 0.5*( bc*bc*fba - ba*ba*fbc )/( bc*fba - ba*fbc );
        if ( std::abs((u-b)/b) < frac_tol ) likelyConvergedToMinimum = true;
        return;
      }
      else if ( fu >= fb ) {
        //bracket found: (a,b,c) <- (a,b,u)
        c = u; fc = fu;
        foundBracket = true;
        //cout << " a " << a << " b " << b << " c " << c <<
          //" fa " << fa << " fb " << fb << " fc " << fc << endl;
        return;
      }
      else { // fu lies between fb and fc
        //u is of no help in
        //setting the next triple.
        //Set next triple by magnifying c.
        a = b; fa = fb;
        b = c; fb = fc;
        c = b + beta*(b-a);
        //cout << "Lambda c (magnification of a,b) " << c << endl;
        fc = f(c);
      }
    }
    else if ( ( b<c && u>c && u<=limit ) ||
              ( b>c && u<c && u>=limit ) ) {
      //We can use u in the next triple
      //Set next triple to (a,b,c) <- (b,c,u)
      //cout << "Lambda u (u>c) " << u << endl;
      fu = f(u);
      a = b; fa = fb;
      b = c; fb = fc;
      c = u; fc = fu;
    }
    else if ( ( b<c && u>c && u>limit ) ||
              ( b>c && u<c && u<limit ) ) {
      //u exceeds the limit.
      //Set next triple to (a,b,c) <- (b,c,limit)
      a = b; fa = fb;
      b = c; fb = fc;
      c = limit;
      //cout << "Lambda c (c=limit) " << c << endl;
      fc = f(limit);
    }
    else { // can happen if u<b, or if u==b, or if u==c
      //u is of no help in
      //setting the next triple.
      //Set next triple by magnifying c.
      a = b; fa = fb;
      b = c; fb = fc;
      c = b + beta*(b-a);
      //cout << "Lambda c (magnifying) " << c << endl;
      fc = f(c);
    }
  }
}

#endif /*FINDBRACKET_H_*/
