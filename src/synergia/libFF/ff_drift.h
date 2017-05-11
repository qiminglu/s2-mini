#ifndef FF_DRIFT_H
#define FF_DRIFT_H

#include "synergia/libFF/ff_algorithm.h"
#include "synergia/libFF/ff_element.h"

class FF_drift : public FF_element
{
private:
    double get_reference_cdt(double length, Reference_particle & reference_particle);
public:
    FF_drift();
    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);
    virtual ~FF_drift();
};

typedef std::shared_ptr<FF_drift > FF_drift_sptr;

#endif // FF_DRIFT_H
