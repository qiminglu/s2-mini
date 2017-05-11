#ifndef BUNCH_SIMULATOR_H_
#define BUNCH_SIMULATOR_H_

#include "synergia/bunch/bunch.h"

class Simulator
{
private:
    Bunch & bunch;

public:
    Simulator(Bunch & bunch) : bunch(bunch) { }
    Bunch & get_bunch() { return bunch; }
};
#endif /* BUNCH_SIMULATOR_H_ */
