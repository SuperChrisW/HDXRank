#ifndef _PHASESPACEPATHLENGTHINTEGRATOR_H_
#define _PHASESPACEPATHLENGTHINTEGRATOR_H_
#include <vector>
#include <cmath>

class PerturbRelaxCycle;
#include "Observable.h"
#include "Vec3.h"

////////////////////////////////////////////////////////////////////////////////
//
//   Scott Menor
//   Arizona State University
//   Department of Physics and Astronomy
//   Biophysics Theory Group
//   PO Box 871504
//   Tempe, AZ 85287
//   Scott.Menor@asu.edu
//
//   Copyright (c) 2005-2007, Arizona State University, All Rights Reserved.
////////////////////////////////////////////////////////////////////////////////
class PhaseSpacePathLengthIntegrator : public Observer {
public:
  PhaseSpacePathLengthIntegrator( std::string pdbfilename,
    PerturbRelaxCycle *cycle_,
    const std::vector<Vec3> *coordinates_,
    double timeLimitPS_ );
  ~PhaseSpacePathLengthIntegrator();

  void activateMassWeighting(const std::vector<double>& masses_ );
  void activateMasking( const std::vector<char>& boolMask_ );

  void updatePath();

  void setTimePerDistance( double timePerDistance_ ) { timePerDistance = timePerDistance_; }
  double getIntegratedPathLength() const { return integratedPathLength; }
  double getIntegratedPathTime() const { return integratedPathLength*timePerDistance; }
  double getSegmentLength() const { return length; }
  void receiveNotification( Observable *obs );


private:
  const std::vector<Vec3> *coordinates;
  std::vector<Vec3> previousCoordinates;
  std::vector<double> masses;
  std::vector<char> boolMask;
  double length;
  double integratedPathLength;
  bool weighting;
  bool masking;
  double timePerDistance;
  double timeLimitPS;

  double computeDistance(
      const std::vector<Vec3> &initialCoordinates,
      const std::vector<Vec3> &finalCoordinates );
};

#endif // _PHASESPACEPATHLENGTHINTEGRATOR_H_
