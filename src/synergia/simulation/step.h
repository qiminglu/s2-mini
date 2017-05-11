#ifndef STEP_H_
#define STEP_H_

#include <list>
#include "synergia/utils/logger.h"
#include "synergia/simulation/operator.h"
#include "synergia/bunch/bunch.h"


class Step
{
private:
    std::vector<std::shared_ptr<Operator>> operators;
    std::list<double> time_fractions;
    double length;
   

public:
    Step(double length) : operators(), time_fractions(), length(length) { }

#if 0
    void append(Operator const & op, double time_fraction)
    { operators.push_back(op); time_fractions.push_back(time_fraction); }
#endif

    Operator * append(std::shared_ptr<Operator> op, double time_fraction);

    void apply(Bunch & bunch, int verbosity, Logger & logger);
    
    std::vector<std::shared_ptr<Operator>> const & get_operators() const { return operators; }
    std::vector<std::shared_ptr<Operator>>       & get_operators()       { return operators; }

    std::list<double>     const & get_time_fractions() const { return time_fractions; }
    double                        get_length()         const { return length; }
};

#endif /* STEP_H_ */
