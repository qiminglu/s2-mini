#include <iostream>
#include <stdexcept>

#include "synergia/lattice/lattice.h"
#include "synergia/simulation/operator.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/propagator.h"
#include "synergia/bunch/bunch.h"
#include "synergia/foundation/physical_constants.h"


// We put the actual code in a separate function so that shared_ptr's can
// be cleanup up properly before we call MPI_Finalize.
void run()
{
    std::vector<int > grid_shape = { 32, 32, 32 };

    const int part_per_cell = 10;
    const int num_macro_particles = grid_shape[0] * grid_shape[1] * grid_shape[2] * part_per_cell;

    const int seed = 4;
    const double num_real_particles = 1e13;

    const int num_steps = 2;
    const int num_turns = 4;
    const int map_order = 2;

    int verbosity = 2;

#if 0
    Lattice_sptr lattice_sptr(new Lattice());
    lattice_sptr->set_all_string_attribute("extractor_type","chef_propagate");
#endif

    Commxx_sptr comm = std::make_shared<Commxx>();

    // reference particle
    Reference_particle ref_part(1, pconstants::mp, 7.0);

    // lattice elements
    Lattice_element e_drift("drift", "d1");
    e_drift.set_double_attribute("l", 1.0);

    // the lattice
    Lattice lattice;
    lattice.set_reference_particle(ref_part);
    lattice.append_element(e_drift);

    Bunch bunch(ref_part, num_macro_particles, num_real_particles, comm);

#if 0
    Random_distribution distribution(seed, *comm_sptr);
    MArray1d means;
    xml_load(means, "cxx_means.xml");
    MArray2d covariances;
    xml_load(covariances, "cxx_covariance_matrix.xml");
    populate_6d(distribution, *bunch_sptr, means, covariances);
    if (opts.sortperiod > 0) {
        bunch_sptr->sort(Bunch::z);
    }
#endif

#if 0
    Space_charge_3d_open_hockney_sptr space_charge_sptr(new Space_charge_3d_open_hockney(grid_shape));
#endif

    //Space_charge_3d_open_hockney space_charge();
    Dummy_collective_operator space_charge("dummy");
    Stepper stepper(lattice, map_order, space_charge, num_steps);

#if 0
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    Split_operator_stepper_sptr stepper_sptr(
            new Split_operator_stepper(lattice_simulator, space_charge_sptr,
                    num_steps));
#endif

    Propagator propagator(stepper);
    Simulator simulator(bunch);

    double t0 = MPI_Wtime();

    const int max_turns = 0;
    propagator.propagate(simulator, num_turns, max_turns, verbosity);

    double t1 = MPI_Wtime();
    std::cout << "propagate time = " << (t1 - t0) << std::endl;
}

int
main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    run();

    MPI_Finalize();
    return 0;
}