
#ifndef LATTICE
#define LATTICE

#include "lattice_element.h"
#include "synergia/foundation/reference_particle.h"


class Lattice
{
public:

    Reference_particle const & get_reference_particle() const { return ref_part; }
    Reference_particle       & get_reference_particle()       { return ref_part; }

    std::vector<Lattice_element> & get_elements() { return vle; }

private:

    std::vector<Lattice_element> vle;
    Reference_particle ref_part;
};


#endif
