#ifndef STEPPER_H_
#define STEPPER_H_

#include <list>

#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/operator.h"
#include "synergia/simulation/step.h"

class Stepper
{
public:
    static const std::string force_diagnostics_attribute;
    static const double fixed_step_tolerance;

private:
    Lattice_simulator lattice_simulator;
    std::vector<Step> steps;

protected:
#if 0
    Independent_operator_sptr get_fixed_step(std::string const& name,
        Lattice_elements::iterator & lattice_it, double & left,
        Lattice_elements::iterator const & lattice_end,
        const double step_length, double & offset_fudge,
        bool end_on_force_diagnostics);

    Lattice_element_slices extract_slices(Steps const& steps);
#endif

public:
    Stepper(Lattice & lattice, 
            int map_order, 
            Collective_operator & col_op, 
            int steps_per_element);

    Lattice_simulator & get_lattice_simulator()
    { return lattice_simulator; }

    std::vector<Step> & get_steps()
    { return steps; }

#if 0
    void force_update_operations_no_collective();
#endif
};

#endif /* STEPPER_H_ */
