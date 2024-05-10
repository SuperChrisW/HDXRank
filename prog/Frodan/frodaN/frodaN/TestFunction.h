#ifndef TESTFUNCTION_H_
#define TESTFUNCTION_H_

class TestFunction
{
public:
	TestFunction() {}
	virtual ~TestFunction() {}
  double operator()( double x ) {
    return (x-5.0)*(x-5.0);
  }
};

#endif /*TESTFUNCTION_H_*/
