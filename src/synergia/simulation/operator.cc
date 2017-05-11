#include <iostream>
#include "synergia/utils/simple_timer.h"
#include "operator.h"

void Independent_operator::update_operations( Reference_particle const& reference_particle)
{
    operations.clear();
    operations_revisions.clear();

    std::string aperture_type("");
    std::string extractor_type(""), last_extractor_type("");

    // Group slices of equal extractor_type and pass to operation_extractor
    // to get operations.
    extractor_type = "libff";
    operations.push_back(std::unique_ptr<Independent_operation>(new LibFF_operation(slices)));

#if 0
    Lattice_element_slices group;
    for (Lattice_element_slices::iterator it = slices.begin(); it != slices.end(); ++it) 
    {
        if ((*it)->get_lattice_element().has_string_attribute("aperture_type")) 
        {
            aperture_type = (*it)->get_lattice_element().get_string_attribute("aperture_type");
            need_left_aperture = (*it)->has_left_edge();
            need_right_aperture = (*it)->has_right_edge();
        } 
        else 
        {
            need_left_aperture = false;
            need_right_aperture = false;
        }

        if ((*it)->get_lattice_element().has_string_attribute("extractor_type")) 
        {
            extractor_type = (*it)->get_lattice_element().get_string_attribute("extractor_type");
        } 
        else 
        {
            extractor_type = "default";
        }

        if (((extractor_type != last_extractor_type) || need_left_aperture) && (!group.empty())) 
        {
            Independent_operations group_operations =
                    operation_extractor_map_sptr->get_extractor(last_extractor_type)->extract(
                            reference_particle, group);
            operations.splice(operations.end(), group_operations);
            group.clear();
        }

        if (need_left_aperture) 
        {
            Aperture_operation_extractor_sptr extractor(
                    aperture_operation_extractor_map_sptr->get_extractor(aperture_type));
            Aperture_operation_sptr aperture_operation_sptr(
                    extractor->extract(*it));
            operations.push_back(aperture_operation_sptr);
        }

        group.push_back(*it);
        last_extractor_type = extractor_type;

        if (need_right_aperture) 
        {
            Aperture_operation_extractor_sptr extractor(
                    aperture_operation_extractor_map_sptr->get_extractor(
                            aperture_type));
            Aperture_operation_sptr aperture_operation_sptr(
                    extractor->extract(*it));
            operations.push_back(aperture_operation_sptr);
            Independent_operations group_operations =
                    operation_extractor_map_sptr->get_extractor(extractor_type)->extract(
                            reference_particle, group);
            operations.splice(operations.end(), group_operations);
            group.clear();
        }

        operations_revisions.push_back((*it)->get_lattice_element().get_revision());
    }

    if (!group.empty()) 
    {
        Independent_operations group_operations =
                operation_extractor_map_sptr->get_extractor(extractor_type)->extract(reference_particle, group);
        operations.splice(operations.end(), group_operations);
    }

    Aperture_operation_extractor_sptr extractor(aperture_operation_extractor_map_sptr->get_extractor("default"));
    Aperture_operation_sptr aperture_operation_sptr(extractor->extract(slices.back()));
    operations.push_back(aperture_operation_sptr);

    Aperture_operation_sptr finite_aperture_operation_sptr(new Finite_aperture_operation(slices.back()));
    operations.push_back(finite_aperture_operation_sptr);
#endif

    have_operations = true;
    operations_reference_particle = reference_particle;
}

bool Independent_operator::need_update(Reference_particle const& reference_particle, int verbosity, Logger & logger)
{
    const double reference_particle_tolerance = 1.0e-8;
    bool retval;

    if (have_operations) 
    {
        retval = false;

        if (reference_particle.equal(operations_reference_particle, reference_particle_tolerance)) 
        {
            retval = false;
#if 0
            std::list<long int>::const_iterator rev_it = operations_revisions.begin();
            for (Lattice_element_slices::const_iterator it = slices.begin(); it != slices.end(); ++it) 
            {
                long int cached_revision = (*rev_it);
                long int revision = (*it)->get_lattice_element().get_revision();

                if (revision != cached_revision) 
                {
                    if (verbosity > 4) 
                    {
                        logger << "Independent_operator: needs update because lattice element "
                               << (*it)->get_lattice_element().get_name()
                               << " has changed" << std::endl;
                    }

                    retval = true;
                    break;
                }

                ++rev_it;
            }
#endif
        } 
        else 
        {
            if (verbosity > 4) logger << "Independent_operator: needs update because reference particle has changed" << std::endl;
            retval = true;
        }
    } 
    else 
    {
        if (verbosity > 4) logger << "Independent_operator: needs update because does not have operations" << std::endl;
        retval = true;
    }

    return retval;
}

void
Independent_operator::apply(Bunch & bunch, double time_step, Step & step, int verbosity, Logger & logger)
{
    double t_total = simple_timer_current();
    double t = simple_timer_current();

    bool do_update = need_update(bunch.get_reference_particle(), verbosity, logger);
    t = simple_timer_show(t, "independent_operator_apply-test_update");

    if (do_update) 
    {
        if (verbosity > 3) logger << "Independent_operator: updating operations" << std::endl;

        update_operations(bunch.get_reference_particle());
        t = simple_timer_show(t, "independent_operator_apply-update_operations");
    }

    for (auto const & opn : operations)
    {
        if (verbosity > 3)
            logger << "Independent_operator: operation type = " << opn->get_type() << std::endl;

        opn->apply(bunch, verbosity, logger);

        std::string label("independent_operator_apply-" + opn->get_type() + "_operation_apply");
        t = simple_timer_show(t, label.c_str());
    }

    bunch.update_total_num();

    t = simple_timer_show(t, "independent_operator_apply-bunch_update_total_num");
    t_total = simple_timer_show(t_total, "independent_operator_apply-total");
}

