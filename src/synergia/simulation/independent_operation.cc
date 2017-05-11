#include "independent_operation.h"
#include "synergia/libFF/ff_element_map.h"
#include "synergia/utils/simple_timer.h"


void LibFF_operation::apply(Bunch & bunch, int verbosity, Logger & logger)
{
    double t = simple_timer_current();

    bunch.convert_to_state(Bunch::fixed_z_lab);
    t = simple_timer_show(t, "LibFF_operation_apply-convert_to_state");

    for (auto const & slice : lattice_element_slices)
    {
        std::string const & type(slice->get_lattice_element().get_type());
        the_big_giant_global_ff_element_map.get_element_type(type)->apply(*slice, bunch);
    }

    t = simple_timer_show(t, "LibFF_operation_apply-element_apply");
}
