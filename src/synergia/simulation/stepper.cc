#include "stepper.h"
#include "synergia/utils/floating_point.h"
//#include "synergia/utils/string_utils.h"
#include "synergia/bunch/bunch.h"
#include <cmath>

#if 0
const std::string Stepper::force_diagnostics_attribute("force_diagnostics");
const double Stepper::fixed_step_tolerance = 1.0e-8;
#endif

Stepper::Stepper(Lattice & lattice, int map_order, Collective_operator & col_op, int steps_per_element) 
: lattice_simulator(lattice, map_order)
, steps()
{
    if (steps_per_element < 1) 
    {
        throw std::runtime_error("Split_operator_stepper_elements: steps_per_element must be >= 1");
    }

    // lattice elements
    auto & lattice_elements = lattice_simulator.get_lattice().get_elements();

    // loop through elements
    for (auto & le : lattice_elements)
    {
        le.clear_slices();
        double length = le.get_length();

        //zero-length element
        if (length == 0.0) 
        {
            Step step(0.0);

            auto ind_op_ptr = std::make_shared<Independent_operator>("step");
            auto slice = le.add_slice();
            dynamic_cast<Independent_operator *>(step.append(ind_op_ptr, 1.0)) -> append_slice(slice);

            steps.push_back(step);
        } 
        else 
        {
            double step_length = length / steps_per_element;
            double half_step_length = 0.5 * step_length;

            for (int i = 0; i < steps_per_element; ++i) 
            {
                double left = i * step_length;
                double middle = left + half_step_length;
                double right = (i + 1) * step_length;

                Step step(step_length);

                //1st Half
                auto ind_op_ptr = std::make_shared<Independent_operator>("first_half");
                auto ind_op_first_half = dynamic_cast<Independent_operator *>(step.append(ind_op_ptr, 0.5));
                auto slice_1st_half = le.add_slice(left, middle);
                ind_op_first_half->append_slice(slice_1st_half);

                //Collective Effects
                //step.append(col_op, 1.0);
#if 0
                for (Collective_operators::const_iterator coll_op_it =
                        collective_operators.begin(); coll_op_it
                        != collective_operators.end(); ++coll_op_it) {
                    Collective_operator_sptr copied_collective_operator_sptr(
                                                            (*coll_op_it)->clone());
                    step->append(copied_collective_operator_sptr, 1.0);
                }
#endif

                //2nd Half
                ind_op_ptr = std::make_shared<Independent_operator>("second_half");
                auto ind_op_second_half = dynamic_cast<Independent_operator *>(step.append(ind_op_ptr, 0.5));
                auto slice_2nd_half = le.add_slice(middle, right);
                ind_op_second_half->append_slice(slice_2nd_half);

                // push to steps
                steps.push_back(step);

            }
        }
    }

    //lattice_simulator.set_slices( extract_slices(steps) );
}

#if 0
void
Stepper::force_update_operations_no_collective()
{
    int total_num = 1;
    double real_num = 1.0;
    Commxx_sptr commxx_sptr(new Commxx);
    Bunch bunch(lattice_simulator.get_lattice_sptr()->get_reference_particle(),
            total_num, real_num, commxx_sptr);
    int verbosity = 0;
    Logger logger(0);
    double time_step = 1.0;
    for (int i = 0; i < 6; ++i) {
        bunch.get_local_particles()[0][i] = 0.0;
    }
    Diagnosticss dummy_diagnosticss;
    for (Steps::iterator sit = steps.begin(); sit != steps.end(); ++sit) {
        for (Operators::iterator oit = (*sit)->get_operators().begin();
                oit != (*sit)->get_operators().end(); ++oit) {
            if ((*oit)->get_type() == Independent_operator::type_name) {
                (*oit)->apply(bunch, time_step, **sit, verbosity,
                        dummy_diagnosticss, logger);
            }
        }
    }
}
#endif

#if 0
// Return an Independent_operator for a half step, starting at the
// lattice_element given by lattice_it at position left. Both lattice_it
// and left are updated by the function.
Independent_operator_sptr
Stepper::get_fixed_step(std::string const& name,
        Lattice_elements::iterator & lattice_it, double & left,
        Lattice_elements::iterator const & lattice_end,
        const double step_length, double & offset_fudge,
        bool end_on_force_diagnostics)
{
    Independent_operator_sptr retval(
            new Independent_operator(name,
                    get_lattice_simulator().get_operation_extractor_map_sptr(),
                    get_lattice_simulator().get_aperture_operation_extractor_map_sptr()));
    double length = offset_fudge;
    bool complete = false;
    while (!complete) {
        double right = (*lattice_it)->get_length();
        if (floating_point_leq(length + (right - left), step_length,
                fixed_step_tolerance)) {
            // The rest of the element fits in the half step
            Lattice_element_slice_sptr slice(
                    new Lattice_element_slice(*lattice_it, left, right));
            retval->append_slice(slice);
            length += (right - left);
            if (end_on_force_diagnostics && slice->has_right_edge()
                    && (*lattice_it)->has_string_attribute(
                            force_diagnostics_attribute)) {
                if (!false_string(
                        (*lattice_it)->get_string_attribute(
                                force_diagnostics_attribute))) {
                    complete = true;
                }
            }
            ++lattice_it;
            left = 0.0;
            if (floating_point_equal(length, step_length,
                    fixed_step_tolerance)) {
                if ((lattice_it == lattice_end)
                        || ((*lattice_it)->get_length() != 0.0)) {
                    complete = true;
                }
            } else {
                if (lattice_it == lattice_end) {
                    throw(std::runtime_error(
                            "get_step stepped beyond end of lattice"));
                }
            }
        } else {
            // Need to take a portion of the element...
            bool end_within_error = false;
            double old_right = right;
            right = step_length - length + left;
            if ((old_right - right) < fixed_step_tolerance) {
                // ... unless we are within an accumulated tolerance of the end
                right = old_right;
                end_within_error = true;
            }
            Lattice_element_slice_sptr slice(
                    new Lattice_element_slice(*lattice_it, left, right));
            retval->append_slice(slice);
            length += (right - left);
            if (end_within_error) {
                ++lattice_it;
                left = 0.0;
            } else {
                left = right;
            }
            complete = true;
        }
    }
    offset_fudge = length - step_length;
    return retval;
}
#endif


#if 0
Lattice_element_slices
Stepper::extract_slices(Steps const& steps)
{
    Lattice_element_slices all_slices;

    for (auto const & step : steps)
    {
        auto const & step_operators = step.get_operators();
        for (auto const & op : step_operators)
        {
            if (op.get_type() == "independent")
            {
            }

#if 0
            if ((*o_it)->get_type() == "independent") {
                Lattice_element_slices
                        element_slices(
                                boost::static_pointer_cast<Independent_operator >(
                                        *o_it)->get_slices());
                all_slices.splice(all_slices.end(), element_slices);
            }
#endif
        }
    }

    return all_slices;
}
#endif




