#include <iostream>
#include "step.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/foundation/physical_constants.h"


void Step::apply(Bunch & bunch, int verbosity, Logger & logger)
{
    double t_total = simple_timer_current();

    std::list<double >::const_iterator fractions_it = time_fractions.begin();
    for (auto const & op : operators)
    {
        // time [s] in accelerator frame
        double time = length / (bunch.get_reference_particle().get_beta() * pconstants::c);        
        double t0 = MPI_Wtime();
        double t = simple_timer_current();

        op->apply(bunch, (*fractions_it) * time, *this, verbosity, logger);

        std::string label("step_apply-" + op->get_type() + "_operator_apply");
        t = simple_timer_show(t, label.c_str());

        double t1 = MPI_Wtime();

        if (verbosity > 2) 
        {
            logger << "Step: operator: name = " << op->get_name()
                    << ", type = " << op->get_type() << ", time = "
                    << std::fixed << std::setprecision(3) << t1 - t0 << "s_n"
                    << std::endl;
        }

        t = simple_timer_current();
        t = simple_timer_show(t, "diagnostics-operator");

        ++fractions_it;
    }

    t_total = simple_timer_show(t_total, "step_apply-total");
}


Operator * Step::append(std::shared_ptr<Operator> op, double time_fraction)
{ 
    operators.push_back(op);
    time_fractions.push_back(time_fraction);
    return operators.back().get();
}


