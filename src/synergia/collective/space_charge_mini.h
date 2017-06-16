#ifndef SPACE_CHARGE_MINI_H
#define SPACE_CHARGE_MINI_H

#include "synergia/simulation/operator.h"
#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/foundation/reference_particle.h"
#include "synergia/utils/commxx.h"

class Space_charge_mini : public Collective_operator
{
public:
    Space_charge_mini(std::string const& name) : Collective_operator(name) 
    { }

    virtual ~Space_charge_mini() 
    { }

    virtual void apply(Bunch & bunch, double time_step, Step & step, int verbosity, Logger & logger);
};

#endif
