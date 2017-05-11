
#ifndef LATTICE_SIMULATOR
#define LATTICE_SIMULATOR

#include "synergia/lattice/lattice.h"

class Lattice_simulator
{
public:

    Lattice_simulator(Lattice & lattice, int map_order)
      : lattice(lattice), map_order(map_order)
    { }

    Lattice & get_lattice() { return lattice; }

private:

    Lattice & lattice;
    int map_order;

};


#endif
