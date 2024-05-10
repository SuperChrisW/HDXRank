#ifndef LINEMINIMIZER_H_
#define LINEMINIMIZER_H_

#include <cmath>
using namespace std;

template <typename T>
class LineMinimizer
{
public:
	LineMinimizer();
	virtual ~LineMinimizer() {}
  double operator()( T &f, double tol, 
                    double a, double b, double c,
                    double fa, double fb, double fc );
  int count;
private:
  const double R1;
  const double R2;
};

template <typename T>
LineMinimizer<T>::LineMinimizer() :
  R1(0.5*(3.0 - sqrt(5.0))),
  R2(1-R1) {}

template <typename T>
double LineMinimizer<T>::operator()( T &f, double tol, 
      double a, double b, double c,
      double fa, double fb, double fc )
{
  double d;
  double temp;
  
  // If the input bracket is already tight enough,
  // we are done.
  if ( (c-a) <= tol ) return b;
  
  // create ordered set a,b,c,d with the next trial value at either b or c
  if ( (b-a) < (c-b) ) {
    //put original bracket a,b,c into variables a,b,d
    //and set c to the trial value between b and d
    d = c;
    c = b + R1*(d-b); fc = f(c);
  }
  else {
    //put original bracket a,b,c into variables a,c,d
    //and set b to the trial value between a and c
    d = c;
    c = b; fc = fb;
    b = c - R1*(c-a); fb = f(b);
  }  

  /*
  if ( (b-a) < (c-b) ) {
    a = bracketa;
    b = bracketb;
    c = bracketb + R1*(bracketc-bracketb);
    d = bracketc;
  }
  else {
    a = bracketa;
    b = bracketb - R1*(bracketb-bracketa);
    c = bracketb;
    d = bracketc;
  }  
  fb = f(b);
  fc = f(c);
  */
  count = 0;
  //shrink the bracket until the tolerance condition is met
  while ( (d-a) > tol ) {
    count++;
    // If f(c) > f(b), the new bracket with new trial point is 
    //    (a',b',c',d') <- (a,b-R1*(b-a),b,c)
    // Otherwise, the new bracket with new trial point is
    //    (a',b',c',d') <- (b,c,c+R1*(d-c),d)
    if ( fc > fb ) {
      temp = b-R1*(b-a);
      d = c;
      c = b; fc = fb;
      b = temp; fb = f(b);
    }
    else {
      temp = c+R1*(d-c);
      a = b;
      b = c; fb = fc;
      c = temp; fc = f(c);
    }
  }
  return ( fb < fc ) ? b : c;  
}

#endif /*LINEMINIMIZER_H_*/
