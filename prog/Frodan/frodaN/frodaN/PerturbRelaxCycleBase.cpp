#include "PerturbRelaxCycleBase.h"

using namespace std;

PerturbRelaxCycleBase::PerturbRelaxCycleBase() :
  count(0),
  state(0),
  go( true ),
  stopAtNsteps(0)
{
}

PerturbRelaxCycleBase::~PerturbRelaxCycleBase()
{
}

void PerturbRelaxCycleBase::doCycle()
{
  count++;
  state = 0;
  notifyObservers(); //state is 0

  commands_cycleStart();
  state++;
  notifyObservers(); //state is 1

  commands_perturb();
  state++;
  notifyObservers(); //state is 2

  commands_minimize();
  state++;
  notifyObservers(); //state is 3

  commands_postMinimization();
  state++;
  notifyObservers(); //state is 4
}
