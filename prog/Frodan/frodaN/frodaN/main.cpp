#include "Settings.h"
#include "PerturbRelaxCycle.h"
#include "Timer.h"
#include "mt19937ar.h"
#include "RigidUnitSystem.h"
#include "Targeter.h"
#include "DynamicConstraints.h"
#include <iostream>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

int main( int argc, char **argv ) {
  #ifdef _OPENMP
  cout << "Compiled with OpenMP; " << omp_get_num_procs() << " processors available.\n";
  #endif

  Settings settings( argc, argv );

  cout << hex << showbase << internal << setfill('0');
  cout << "Seed for Random Number Generator: " << setw(10) << settings.seed << endl;
  cout << dec << noshowbase;
  init_genrand(settings.seed); //declared in mt19937ar.h

  if ( settings.runtype == "target" ) {
    Targeter* sim = new Targeter( settings );
    sim->run();
    delete sim;
  }
  else if ( settings.runtype == "fixedcon" || settings.runtype == "dyncon") {
    //Note, even though this object is called "Dynamic Constraints", it
    //has the capability of running in either "fixed constraints" or
    //"dynamic constraints" mode.
    DynamicConstraints* sim = new DynamicConstraints( settings );
    sim->run( settings.perturbRelax.Nsteps );
    delete sim;
  }
  else if ( settings.runtype == "fixedcon_old" ) {
    PerturbRelaxCycle* sim = new PerturbRelaxCycle( settings );
    sim->doCycles();
    delete sim;
  }
  else {
    cout << "No valid run type specified." << endl;
  }

  return(0);
}
