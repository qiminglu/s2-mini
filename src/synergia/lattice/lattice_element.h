#ifndef LATTICE_ELEMENT_H_
#define LATTICE_ELEMENT_H_

#include <string>
#include <map>
#include <list>
#include <vector>

class Lattice;
class Lattice_element;

class Lattice_element_slice
{
private:
    Lattice_element * element_ptr;
    bool whole;
    bool left_edge;
    bool right_edge;
    double left;
    double right;

public:
    Lattice_element_slice(Lattice_element & element);
    Lattice_element_slice(Lattice_element & element, double left, double right);

    bool   is_whole()       const { return whole; }
    bool   has_left_edge()  const { return left_edge; }
    bool   has_right_edge() const { return right_edge; }
    double get_left()       const { return left; }
    double get_right()      const { return right; }

    Lattice_element const & get_lattice_element() const { return *element_ptr; }
    Lattice_element       & get_lattice_element()       { return *element_ptr; }
};



/// The Lattice_element class contains the description of a single
/// lattice element. Each element has a name, a (string) type and
/// dictionaries of named double and string attributes.
/// Lattice structure is described by a list of ancestors stored in
/// an element.
class Lattice_element
{
private:
    std::string type;
    std::string name;

    std::map<std::string, double>              double_attributes;
    std::map<std::string, std::string>         string_attributes;
    std::map<std::string, std::vector<double>> vector_attributes;

    std::string length_attribute_name;
    std::string bend_angle_attribute_name;

    long int revision;
    bool needs_internal_derive, needs_external_derive;

    Lattice * lattice_ptr;

    std::vector<Lattice_element_slice> slices;

public:
    /// Construct a Lattice_element with an empty name and type.
    Lattice_element();

    /// Construct a Lattice_element.
    /// @param name name
    /// @param type type
    Lattice_element(std::string const& type, std::string const& name);

    /// Get the type
    std::string const & get_type() const { return type; }

    /// Get the name
    std::string const & get_name() const { return name; }

    /// addons
    void clear_slices() 
    { slices.clear(); }

    Lattice_element_slice & add_slice() 
    { slices.emplace_back(*this); return slices.back(); }

    Lattice_element_slice & add_slice(double left, double right)
    { slices.emplace_back(*this, left, right); return slices.back(); }

    /// Set the value of the named double attribute
    /// @param name attribute name
    /// @param value attribute value
    /// @param increment_revision can be set to false for attributes that do not affect dynamics
    void
    set_double_attribute(std::string const& name, double value, bool increment_revision = true);

    /// Check for the existence of the named double attribute
    /// @param name attribute name
    bool
    has_double_attribute(std::string const& name, bool include_default = true) const;

    /// Get the value of the named double attribute
    /// @param name attribute name
    double
    get_double_attribute(std::string const& name) const;

    /// Get the value of the named double attribute
    /// @param name attribute name
    /// @param val default value if the specified attribute doesnt exist
    double
    get_double_attribute(std::string const& name, double val) const;

    /// Get the entire dictionary of double attributes
    std::map<std::string, double > const &
    get_double_attributes() const;

    /// Set the value of the named string attribute
    /// @param name attribute name
    /// @param value attribute value
    /// @param increment_revision can be set to false for attributes that do not affect dynamics
    void
    set_string_attribute(std::string const& name, std::string const& value, bool increment_revision = true);

    /// Check for the existence of the named string attribute
    /// @param name attribute name
    bool
    has_string_attribute(std::string const& name, bool include_default = true) const;

    /// Get the value of the named string attribute
    /// @param name attribute name
    std::string const&
    get_string_attribute(std::string const& name) const;

    /// Get the value of the named string attribute
    /// @param name attribute name
    /// @param val default value if the specified attribute doesnt exist
    std::string const&
    get_string_attribute(std::string const& name, std::string const & val) const;

    /// Get the entire dictionary of string attributes
    std::map<std::string, std::string > const &
    get_string_attributes() const;

    /// Set the value of the named vector attribute
    /// @param name attribute name
    /// @param value attribute value
    /// @param increment_revision can be set to false for attributes that do not affect dynamics
    void
    set_vector_attribute(std::string const& name, std::vector<double > const& value, bool increment_revision = true);

    /// Check for the existence of the named vector attribute
    /// @param name attribute name
    bool
    has_vector_attribute(std::string const& name, bool include_default = true) const;

    /// Get the value of the named vector attribute
    /// @param name attribute name
    std::vector<double > const&
    get_vector_attribute(std::string const& name) const;

    /// Get the value of the named vector attribute
    /// @param name attribute name
    /// @param val default value if the specified attribute doesnt exist
    std::vector<double > const&
    get_vector_attribute(std::string const& name, std::vector<double> const & val) const;

    /// Get the entire dictionary of vector attributes
    std::map<std::string, std::vector<double > > const &
    get_vector_attributes() const;

    /// Set the attribute name to be used to determine the length
    /// of the Lattice_element
    /// @param attribute_name attribute name
    void
    set_length_attribute_name(std::string const& attribute_name);

    /// Set the attribute name to be used to determine the bend_angle
    /// of the Lattice_element
    /// @param attribute_name attribute name
    void
    set_bend_angle_attribute_name(std::string const& attribute_name);

    /// Set whether the element needs to determine some of its parameters
    /// from its other parameters
    void
    set_needs_internal_derive(bool value);

    /// Get whether the element needs to determine some of its parameters
    /// from its other parameters
    bool
    get_needs_internal_derive() const;

    /// Set whether the element needs to determine some of its parameters
    /// from the lattice length and/or reference particle
    void
    set_needs_external_derive(bool value);

    /// Get whether the element needs to determine some of its parameters
    /// from the lattice length and/or reference particle
    bool
    get_needs_external_derive() const;

    /// Get the Lattice_element's length
    double
    get_length() const;

    /// Get the Lattice_element's bend angle
    double
    get_bend_angle() const;

    /// Get the Lattice_element's revision number
    long int
    get_revision() const;

    /// Check whether the element has a reference to a parent lattice
    bool
    has_lattice() const;

    /// Set the reference to the parent lattice
    void
    set_lattice(Lattice &lattice);

    /// Get a reference to the parent lattice
    Lattice &
    get_lattice();

    /// Get a reference to the parent lattice
    Lattice const&
    get_lattice() const;
};

#endif /* LATTICE_ELEMENT_H_ */
