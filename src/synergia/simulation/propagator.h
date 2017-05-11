#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

#include "synergia/utils/logger.h"
#include "synergia/simulation/stepper.h"
#include "synergia/simulation/simulator.h"

class Propagator
{
public:

    struct State
    {
        Simulator & simulator;

        int num_turns;
        int first_turn;
        int max_turns;
        int verbosity;

        State(Simulator & simulator, int num_turns, int first_turn, int max_turns, int verbosity)
          : simulator(simulator)
          , num_turns(num_turns)
          , first_turn(first_turn)
          , max_turns(max_turns)
          , verbosity(verbosity)
        { }
    };

private:

    Stepper & stepper;

    int omp_threads;

	void do_before_start(State & state, double & t, Logger & logger);
	void do_step(Step & step, int step_count, int num_steps, int turn, State & state, double & t, Logger & logger);
	void do_start_repetition(State & state);
	void do_turn_end(int turn, State & state, double & t, double t_turn0, Logger & logger);

	bool check_out_of_particles(State & state, Logger & logger);

public:

    Propagator(Stepper & stepper)
      : stepper(stepper)
      , omp_threads(1)
    { }

    ~Propagator()
    { }

    Stepper & get_stepper() 
    { return stepper; }

    void set_num_threads(int nt) 
    { omp_threads = nt; }

    int  get_num_threads() const 
    { return omp_threads; }

    void propagate(State & state);
    void propagate(Simulator & simulator, int num_turns, int max_turns = 0, int verbosity = 1);

};

#endif /* PROPAGATOR_H_ */
