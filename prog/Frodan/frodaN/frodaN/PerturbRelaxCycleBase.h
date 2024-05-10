#ifndef PERTURBRELAXCYCLEBASE_H_
#define PERTURBRELAXCYCLEBASE_H_

#include "Command.h"
#include "CommandList.h"
#include "CommandMemberFunction.h"
#include "Observable.h"

class PerturbRelaxCycleBase : public Observable
{
public:
  PerturbRelaxCycleBase();
  virtual ~PerturbRelaxCycleBase();

  void doCycle( int N ) {
    for ( int i = 0; i < N; i++ ) {
      if ( !go ) break;
      doCycle();
    }
  }

  void doCycles() {
    if ( stopAtNsteps ) {
      for ( int i = 0; i < stopAtNsteps; i++ ) {
        if ( !go ) break;
        doCycle();
      }
    }
    else {
      while ( go ) doCycle();
    }
  }

  void doCycle();
  void stop() { go = false; }
  void setStopAtNSteps( int N ) { stopAtNsteps = N; }
  unsigned int getCycleCount() const { return count; }
  unsigned int getState() const { return state; }

  template <class T> void addCommand_cycleStart( T *t, void (T::*f)() ) {
    commands_cycleStart.addCommand( new CommandMemberFunction<T>( t, f ) );
  }
  template <class T> void addCommand_perturb( T *t, void (T::*f)() ) {
    commands_perturb.addCommand( new CommandMemberFunction<T>( t, f ) );
  }
  template <class T> void addCommand_minimize( T *t, void (T::*f)() ) {
    commands_minimize.addCommand( new CommandMemberFunction<T>( t, f ) );
  }
  template <class T> void addCommand_postMinimization( T *t, void (T::*f)() ) {
    commands_postMinimization.addCommand( new CommandMemberFunction<T>( t, f ) );
  }

private:
  CommandList commands_cycleStart;
  CommandList commands_perturb;
  CommandList commands_minimize;
  CommandList commands_postMinimization;
  unsigned int count;
  unsigned int state;
  bool go;
  int stopAtNsteps;
};

#endif /*PERTURBRELAXCYCLEBASE_H_*/
