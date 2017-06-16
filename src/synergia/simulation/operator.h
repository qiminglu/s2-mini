#ifndef OPERATOR_H_
#define OPERATOR_H_

#include <string>
#include <list>

#include "synergia/bunch/bunch.h"
#include "synergia/lattice/lattice_element.h"
#include "synergia/simulation/independent_operation.h"
#include "synergia/utils/logger.h"


class Step;

// Base Operator
//
class Operator
{
private:
    std::string name, type;

public:
    Operator(std::string const& name, std::string const& type)
      : name(name), type(type)
    { }

    virtual ~Operator() 
    { }

    std::string const& get_name() const { return name; }
    std::string const& get_type() const { return type; }

    virtual void apply(Bunch & bunch, double time_step, Step & step, int verbosity, Logger & logger) = 0;
};

// Collective Operator
//
class Collective_operator : public Operator
{
public:
    Collective_operator(std::string const& name) : Operator(name, "collective") { }
    virtual ~Collective_operator() { }
    virtual void apply(Bunch & bunch, double time_step, Step & step, int verbosity, Logger & logger) { }
};

typedef std::shared_ptr<Collective_operator> Collective_operator_sptr;

// Dummy collective operator
//
class Dummy_collective_operator : public Collective_operator
{
public:
    Dummy_collective_operator(std::string const& name) : Collective_operator(name) { }
    virtual ~Dummy_collective_operator() { }
    virtual void apply(Bunch & bunch, double time_step, Step & step, int verbosity, Logger & logger) { logger<<"dummy apply\n"; }
};


// Independent Operator
//
class Independent_operator : public Operator
{
private:
    std::vector<Lattice_element_slice const *> slices;
    bool have_operations;

    std::vector<std::unique_ptr<Independent_operation>> operations;
    std::list<long int > operations_revisions;
    Reference_particle operations_reference_particle;

public:
    Independent_operator( std::string const& name )
      : Operator(name, "independent")
      , have_operations(false)
      , operations()
      , operations_revisions()
      , operations_reference_particle()
    { }

    virtual ~Independent_operator() { }

    void update_operations(Reference_particle const& reference_particle);
    bool need_update(Reference_particle const& reference_particle, int verbosity, Logger & logger);

    void append_slice(Lattice_element_slice const * slice)
    { slices.push_back(slice); have_operations = false; }

    std::vector<Lattice_element_slice const *> const & get_slices() const { return slices; }

    std::vector<std::unique_ptr<Independent_operation>> const & get_operations() const { return operations; }
    std::vector<std::unique_ptr<Independent_operation>>       & get_operations()       { return operations; }

    virtual void apply(Bunch & bunch, double time_step, Step & step, int verbosity, Logger & logger);
};

#endif /* OPERATOR_H_ */
