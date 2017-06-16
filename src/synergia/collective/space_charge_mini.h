#ifndef SPACE_CHARGE_MINI_H
#define SPACE_CHARGE_MINI_H

#include "synergia/simulation/operator.h"
#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/foundation/reference_particle.h"
#include "synergia/utils/commxx.h"

class Space_charge_mini : public Collective_operator
{
public:
    Space_charge_mini(std::vector<int> const & grid_shape);
    virtual ~Space_charge_mini() { }

    virtual void apply(Bunch & bunch, double time_step, Step & step, int verbosity, Logger & logger);

private:

    void get_charge_density(Bunch & bunch, Logger & logger);

    std::vector<int> shape;
    svec<double> rho;
};

#endif
