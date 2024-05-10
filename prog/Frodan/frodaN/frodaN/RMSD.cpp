#include "RMSD.h"
#include <cmath>

double RMSD::calcRMSD( const std::vector<Vec3>& points1, const std::vector<Vec3>& points2 ) const {
  double dist2sum = 0.0;
  for ( size_t p = 0; p < points1.size(); p++ ) {
    dist2sum += points1[p].dist2( points2[p] );
  }
  return sqrt( dist2sum/static_cast<double>(points1.size()) );
}
