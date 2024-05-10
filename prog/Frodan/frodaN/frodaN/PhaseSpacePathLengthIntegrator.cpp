#include "PhaseSpacePathLengthIntegrator.h"
#include "PerturbRelaxCycle.h"
#include "PDB.h"
#include <cstdlib>

using namespace std;

/**
 *
 */
PhaseSpacePathLengthIntegrator::PhaseSpacePathLengthIntegrator(
    string pdbfilename, PerturbRelaxCycle *cycle, const vector<Vec3> *coordinates_, double timeLimitPS_ ) :
  coordinates( coordinates_ ),
  previousCoordinates( *coordinates ),
  length(0),
  integratedPathLength(0),
  weighting(false),
  masking(false),
  timePerDistance(0),
  timeLimitPS( timeLimitPS_ ) {

  const double timePerPathDistance_psPerAngstrom_CarbonOnly = 1.0/7.894;
  timePerDistance = timePerPathDistance_psPerAngstrom_CarbonOnly;

  PDB pdb( pdbfilename );
  size_t N = pdb.atomLines.size();
  masking = true;
  boolMask.resize( N );
  for ( size_t i = 0; i < N; i++ ) {
    boolMask[i] = pdb.atomLines[i].element == "C";
  }

  cycle->registerObserver( this );
}

/**
 * Destructor
 *
 */
PhaseSpacePathLengthIntegrator::~PhaseSpacePathLengthIntegrator() {
}

void PhaseSpacePathLengthIntegrator::receiveNotification( Observable *obs ) {
  PerturbRelaxCycle *cycle = dynamic_cast<PerturbRelaxCycle*>( obs );
  if ( cycle && cycle->getState() == 4 ) updatePath();
  if ( getIntegratedPathTime() > timeLimitPS ) cycle->stop();
}

void PhaseSpacePathLengthIntegrator::updatePath() {
  length = computeDistance( previousCoordinates, *coordinates );
  integratedPathLength += length;
  previousCoordinates = *coordinates;
}

void PhaseSpacePathLengthIntegrator::activateMassWeighting(const vector<double> &masses_ ) {
  masses = masses_;
  weighting = true;
}

void PhaseSpacePathLengthIntegrator::activateMasking( const vector<char>& boolMask_ ) {
  boolMask = boolMask_;
  masking = true;
}

double PhaseSpacePathLengthIntegrator::computeDistance(
    const vector<Vec3> &initialCoordinates,
    const vector<Vec3> &finalCoordinates ) {

	size_t N = initialCoordinates.size();
	if ( N != finalCoordinates.size() ||
	     ( masking && N != boolMask.size() ) ||
	     ( weighting && N != masses.size() ) ) {
	  cout << "Error in PhaseSpacePathLengthIntegrator: sizes of vectors do not match" << endl;
	  exit(0);
	}
	if (!N) return 0.0;

	// iterate through coordinates and compute
	double sumOfSquaredDistances = 0.0;
	int count = 0;
	double weight = 1.0;
  for (size_t i = 0; i < N; i++) {
    if ( masking && !boolMask[i] ) continue;
    if ( weighting ) weight = masses[i];
    sumOfSquaredDistances += initialCoordinates[i].dist2( finalCoordinates[i] ) * weight;
    count++;
  }

	return sqrt( sumOfSquaredDistances /(double)count );
}
