
#ifndef LATTICE
#define LATTICE

#include "lattice_element.h"
#include "synergia/foundation/reference_particle.h"


class Lattice
{
public:

    void set_reference_particle(Reference_particle const & ref) { ref_part = ref; }

    Reference_particle const & get_reference_particle() const { return ref_part; }
    Reference_particle       & get_reference_particle()       { return ref_part; }

    std::vector<Lattice_element> & get_elements() { return vle; }

    void clear_elements() { vle.clear(); }
    void append_element(Lattice_element const & ele) { vle.push_back(ele); }

private:

    std::vector<Lattice_element> vle;
    Reference_particle ref_part;
};


#endif
