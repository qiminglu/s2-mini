#include <iostream>
#include <stdexcept>
#include <thread>
#include <chrono>

#include "synergia/lattice/lattice.h"
#include "synergia/simulation/operator.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/propagator.h"
#include "synergia/bunch/bunch.h"
#include "synergia/collective/space_charge_mini.h"
#include "synergia/foundation/physical_constants.h"

// We put the actual code in a separate function so that shared_ptr's can
// be cleanup up properly before we call MPI_Finalize.
void run(int ppc)
{
    std::vector<int > grid_shape = { 64, 64, 256 };

    const int part_per_cell = ppc;
    const int num_macro_particles = grid_shape[0] * grid_shape[1] * grid_shape[2] * part_per_cell;

    //const int seed = 4;
    const double num_real_particles = 1e13;

    const int num_steps = 2;
    const int num_turns = 2;
    const int map_order = 2;

    int verbosity = 6;

#if 0
    Lattice_sptr lattice_sptr(new Lattice());
    lattice_sptr->set_all_string_attribute("extractor_type","chef_propagate");
#endif

    Commxx_sptr comm = std::make_shared<Commxx>();

    // reference particle
    Reference_particle ref_part(1, pconstants::mp, 7.0);

    // the lattice
    Lattice lattice;
    lattice.set_reference_particle(ref_part);

    // lattice elements
    Lattice_element e_drift("drift", "d1");
    e_drift.set_lattice(lattice);
    e_drift.set_double_attribute("l", 1.0);

    Lattice_element e_quad("quadrupole", "q1");
    e_quad.set_lattice(lattice);
    e_quad.set_double_attribute("l", 1.0);
    e_quad.set_double_attribute("k1", 0.01);

    // append element to the lattice
    lattice.append_element(e_drift);
    //lattice.append_element(e_quad);

    Bunch bunch(ref_part, num_macro_particles, num_real_particles, comm);
    double * parts = bunch.get_local_particles().origin();
    int npart = bunch.get_local_num();

    for(int i=0; i<npart; ++i)
    {
        parts[npart* 0 + i] = 0.10 + 0.1 * i/num_macro_particles;
        parts[npart* 1 + i] = 0.11 + 0.1 * i/num_macro_particles;
        parts[npart* 2 + i] = 0.12 + 0.1 * i/num_macro_particles;
        parts[npart* 3 + i] = 0.13 + 0.1 * i/num_macro_particles;
        parts[npart* 4 + i] = 0.14 + 0.1 * i/num_macro_particles;
        parts[npart* 5 + i] = 0.15 + 0.1 * i/num_macro_particles;
    }

    for (int i=0; i<7; ++i) std::cout << parts[i*npart + 0] << "\t";
    std::cout << "\n";

    bunch.get_local_particles().copy_host_to_dev();
    //std::this_thread::sleep_for(std::chrono::seconds(5));

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

    auto space_charge = std::make_shared<Space_charge_mini>(grid_shape);
    Stepper stepper(lattice, map_order, space_charge, num_steps);

#if 0
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    Split_operator_stepper_sptr stepper_sptr(
            new Split_operator_stepper(lattice_simulator, space_charge_sptr,
                    num_steps));
#endif

    Simulator simulator(bunch);
    Propagator propagator(stepper);
    propagator.set_num_threads(32);

    double t0 = MPI_Wtime();

    const int max_turns = 0;
    propagator.propagate(simulator, num_turns, max_turns, verbosity);

    double t1 = MPI_Wtime();
    std::cout << "propagate time = " << (t1 - t0) << std::endl;

    bunch.get_local_particles().copy_dev_to_host();

    for (int i=0; i<7; ++i) std::cout << parts[i*num_macro_particles + 0] << "\t";
    std::cout << "\n";
}

int
main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int part_per_cell = 10;

    if (argc == 2)
    {
        part_per_cell = std::stoi(argv[1], nullptr, 0);
    }

    run(part_per_cell);

    MPI_Finalize();
    return 0;
}
