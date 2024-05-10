///////////////////////////////////////////////////////////////////////////////
//
// The source files
//   coeff.h
//   libtricubic.cpp
//   ltricubic_utils.h
//   tricubic.h
//   tricubic_utils.cpp
//
// Are all taken from the tricubic interpolation method of Lekien et al., 
// as described in Int. J. Numer. Meth. Eng. 63:455-471 (2005).
// The code is under a GNU license, and compliled binaries and documentation
// can be found at: http://www.lekien.com/~francois/software/tricubic/
//
///////////////////////////////////////////////////////////////////////////////

char *tricubic_version(void);
void tricubic_get_coeff(double a[64], double f[8], double dfdx[8], double dfdy[8], double dfdz[8], double d2fdxdy[8], double d2fdxdz[8], double d2fdydz[8], double d3fdxdydz[8]);
double tricubic_eval(double a[64], double x, double y, double z);
double tricubic_eval(double a[64], double x, double y, double z, int derx, int dery, int derz);

void tricubic_pointID2xyz(int id, int *x, int *y, int *z);
void tricubic_pointID2xyz(int id, double *x, double *y, double *z);
