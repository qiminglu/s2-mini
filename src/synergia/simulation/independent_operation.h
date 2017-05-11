#ifndef INDEPENDENT_OPERATION_H_
#define INDEPENDENT_OPERATION_H_

#include <list>
#include <string>
#include <map>

#include "synergia/bunch/bunch.h"
#include "synergia/lattice/lattice_element.h"
#include "synergia/utils/logger.h"

class Independent_operation
{
private:
    std::string type;

public:
    Independent_operation(std::string const& type) : type(type) { }
    virtual ~Independent_operation() { }

    std::string const& get_type() const { return type; }

    virtual void apply(Bunch & bunch, int verbosity, Logger & logger) = 0;
};

class LibFF_operation : public Independent_operation
{
private:
    std::vector<Lattice_element_slice const *> const & lattice_element_slices;

public:
    LibFF_operation(std::vector<Lattice_element_slice const *> const & lattice_element_slices)
      : Independent_operation("LibFF"), lattice_element_slices(lattice_element_slices)
    { }

    virtual ~LibFF_operation() 
    { }

    virtual void apply(Bunch & bunch, int verbosity, Logger & logger);
};

#endif /* INDEPENDENT_OPERATION_H_ */
