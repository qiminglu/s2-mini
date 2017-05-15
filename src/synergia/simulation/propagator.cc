
#include "propagator.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/utils/digits.h"

void
Propagator::do_before_start(State & state, double & t, Logger & logger)
{
    if (state.first_turn == 0) 
    {
        //state.simulator.get_diagnostics_actions().first_action(stepper, state.simulator.get_bunch());

        t = simple_timer_show(t, "propagate-diagnostics_actions_first");

        //state.propagate_actions_ptr->first_action(stepper, state.simulator.get_bunch());

        t = simple_timer_show(t, "propagate-general_actions_first");
    }
}

void
Propagator::do_start_repetition(State & state)
{
    state.simulator.get_bunch().get_reference_particle().start_repetition();
}

void
Propagator::do_step(Step & step, int step_count, int num_steps, int turn,
        Propagator::State & state, double & t, Logger & logger)
{
    double t_step0 = MPI_Wtime();

    Bunch & bunch(state.simulator.get_bunch());

    step.apply(bunch, state.verbosity, logger);

    t = simple_timer_show(t, "propagate-step_apply");

    double t_step1 = MPI_Wtime();

    if (state.verbosity > 1) 
    {
        int p = std::cout.precision();
        logger << "Propagator:";
        logger << "     step " << std::setw(digits(num_steps)) << step_count << "/" << num_steps;
        logger << ", s_n = " << std::fixed << std::setprecision(4) << state.simulator.get_bunch().get_reference_particle().get_s_n();
        logger << ", macroparticles = " << state.simulator.get_bunch().get_total_num();
        logger << ", time = " << std::fixed << std::setprecision(3) << t_step1 - t_step0 << "s";
        logger << std::endl;
        std::cout.precision(p);
    }
}

bool
Propagator::check_out_of_particles(State & state, Logger & logger) 
{
    // n.b.: We only check out_of_particles for single-bunch propagation.
    // Checking all bunches in a multi-bunch simulation would require a
    // global communication. The costs exceed the potential benefits.
    if (state.simulator.get_bunch().get_total_num() == 0) 
    {
        logger << "Propagator::propagate: No particles left in bunch. Exiting.\n";
        return true;
    }

	return false;
}

void
Propagator::do_turn_end(int turn, State & state, double & t, double t_turn0, Logger & logger) 
{
    t = simple_timer_current();
    t = simple_timer_show(t, "propagate-general_actions_turn");

    state.first_turn = turn + 1;

    double t_turn1 = MPI_Wtime();

    if (state.verbosity > 0) 
    {
        int p = std::cout.precision();
        logger << "Propagator:";
        logger << " turn " << std::setw(digits(state.num_turns)) << turn + 1 << "/" << state.num_turns;
        logger << ", macroparticles = " << state.simulator.get_bunch().get_total_num();
        logger << ", time = " << std::fixed << std::setprecision(4) << t_turn1 - t_turn0 << "s";
        logger << std::endl;
        std::cout.precision(p);
    }
}

void
Propagator::propagate(State & state)
{
    // set number of openmp threads
    //omp_set_num_threads(omp_threads);

    try 
    {
        Logger logger(0, "logfile");
        double t_total = simple_timer_current();
        double t = simple_timer_current();

        do_before_start(state, t, logger);

        int turns_since_checkpoint = 0;
        bool out_of_particles = false;

        if (state.verbosity > 0) 
            logger << "Propagator: starting turn " << state.first_turn + 1 << std::endl;

        for (int turn = state.first_turn; turn < state.num_turns; ++turn) 
        {
            double t_turn0 = MPI_Wtime();

            do_start_repetition(state);

            int step_count = 0;

            auto steps = stepper.get_steps();
            auto num_steps = steps.size();

            for (auto & step : steps)
            {
                ++step_count;
                do_step(step, step_count, num_steps, turn, state, t, logger);

                out_of_particles = check_out_of_particles(state, logger);
                if (out_of_particles) break;
            }

            if (out_of_particles) break;

            ++turns_since_checkpoint;
            do_turn_end(turn, state, t, t_turn0, logger);	   
        }

        if (out_of_particles) logger << "Propagator: no particles left\n";

        simple_timer_show(t_total, "propagate-total");
    }
    catch (std::exception const& e) 
    {
        std::cerr << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 888);
    }
}

void
Propagator::propagate(Simulator & simulator, int num_turns, int max_turns, int verbosity)
{
    State state(simulator, num_turns, 0, max_turns, verbosity);
    propagate(state);
}


