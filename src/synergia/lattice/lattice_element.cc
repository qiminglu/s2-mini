
#include "lattice_element.h"
#include "lattice.h"
#include "synergia/utils/floating_point.h"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <sstream>

const double split_element_tolerance = 1.0e-9;

Lattice_element_slice::Lattice_element_slice(Lattice_element const & element) 
: element_ptr(&element)
, whole(true), left_edge(true), right_edge(true), left(0.0), right(element.get_length())
{
}

Lattice_element_slice::Lattice_element_slice(Lattice_element const & element, double left, double right) 
: element_ptr(&element), whole(false), left_edge(false), right_edge(false), left(left), right(right)
{
    if (left < 0.0) throw std::range_error("Lattice_element_slice: left must be >= 0.0");

    if (left < split_element_tolerance) 
    {
        left_edge = true;
        left = 0.0;
    } 
    else 
    {
        left_edge = false;
    }

    double element_length = element.get_length();

    if (right > (element_length + split_element_tolerance)) 
        throw std::range_error("Lattice_element_slice: right must be no greater than the length of the element");

    if (floating_point_equal(right, element_length, split_element_tolerance)) 
    {
        right_edge = true;
        right = element_length;
    } 
    else 
    {
        right_edge = false;
    }

    if (left_edge && right_edge) 
    {
        whole = true;
    } 
    else 
    {
        whole = false;
    }
}

Lattice_element::Lattice_element() 
: type("")
, name("")
, double_attributes()
, string_attributes()
, vector_attributes()
, length_attribute_name("l")
, bend_angle_attribute_name("angle")
, revision(0)
, needs_internal_derive(false)
, needs_external_derive(false)
, lattice_ptr(0)
{

}

Lattice_element::Lattice_element(std::string const& type, std::string const& name)
: type("")
, name("")
, double_attributes()
, string_attributes()
, vector_attributes()
, length_attribute_name("l")
, bend_angle_attribute_name("angle")
, revision(0)
, needs_internal_derive(false)
, needs_external_derive(false)
, lattice_ptr(0)
{

}

void
Lattice_element::set_double_attribute(std::string const& name, double value,
        bool increment_revision)
{
    double_attributes[name] = value;
    if (increment_revision) ++revision;
}

bool
Lattice_element::has_double_attribute(std::string const& name,
        bool include_default) const
{
    return double_attributes.find(name) != double_attributes.end();
}

double
Lattice_element::get_double_attribute(std::string const& name) const
{
    std::map<std::string, double >::const_iterator result =
            double_attributes.find(name);
    if (result == double_attributes.end()) {
        throw std::runtime_error("Lattice_element::get_double_attribute: element " + this->name + " of type " + type + " has no double attribute '" + name + "'");
    } else {
        return result->second;
    }
}

double
Lattice_element::get_double_attribute(std::string const& name, double val) const
{
    std::map<std::string, double >::const_iterator result =
            double_attributes.find(name);
    if (result == double_attributes.end()) {
        return val;
    } else {
        return result->second;
    }
}

std::map<std::string, double > const &
Lattice_element::get_double_attributes() const
{
    return double_attributes;
}

void
Lattice_element::set_string_attribute(std::string const& name,
        std::string const& value, bool increment_revision)
{
    string_attributes[name] = value;
    if (increment_revision) {
        ++revision;
    }
}

bool
Lattice_element::has_string_attribute(std::string const& name,
        bool include_default) const
{
    return string_attributes.find(name) != string_attributes.end();
}

std::string const&
Lattice_element::get_string_attribute(std::string const& name) const
{
    std::map<std::string, std::string >::const_iterator result =
            string_attributes.find(name);
    if (result == string_attributes.end()) {
        throw std::runtime_error( "Lattice_element::get_string_attribute: element "
                        + this->name + " of type " + type
                        + " has no string attribute '" + name + "'");
    } else {
        return result->second;
    }
}

std::string const&
Lattice_element::get_string_attribute(std::string const& name, std::string const & val) const
{
    std::map<std::string, std::string >::const_iterator result =
            string_attributes.find(name);
    if (result == string_attributes.end()) {
        return val;
    } else {
        return result->second;
    }
}

std::map<std::string, std::string > const &
Lattice_element::get_string_attributes() const
{
    return string_attributes;
}

void
Lattice_element::set_vector_attribute(std::string const& name,
        std::vector<double > const& value, bool increment_revision)
{
    vector_attributes[name] = value;
    if (increment_revision) {
        ++revision;
    }
}

bool
Lattice_element::has_vector_attribute(std::string const& name,
        bool include_default) const
{
    return vector_attributes.find(name) != vector_attributes.end();
}

std::vector<double > const&
Lattice_element::get_vector_attribute(std::string const& name) const
{
    std::map<std::string, std::vector<double > >::const_iterator result =
            vector_attributes.find(name);
    if (result == vector_attributes.end()) {
        throw std::runtime_error( "Lattice_element::get_vector_attribute: element "
                        + this->name + " of type " + type
                        + " has no vector attribute '" + name + "'");
    } else {
        return result->second;
    }
}

std::vector<double > const&
Lattice_element::get_vector_attribute(std::string const& name, std::vector<double> const & val) const
{
    std::map<std::string, std::vector<double > >::const_iterator result =
            vector_attributes.find(name);
    if (result == vector_attributes.end()) {
        return val;
    } else {
        return result->second;
    }
}

std::map<std::string, std::vector<double > > const &
Lattice_element::get_vector_attributes() const
{
    return vector_attributes;
}

void
Lattice_element::set_length_attribute_name(std::string const& attribute_name)
{
    length_attribute_name = attribute_name;
}

void
Lattice_element::set_bend_angle_attribute_name(
        std::string const& attribute_name)
{
    bend_angle_attribute_name = attribute_name;
}

void
Lattice_element::set_needs_internal_derive(bool value)
{
    needs_internal_derive = value;
}

bool
Lattice_element::get_needs_internal_derive() const
{
    return needs_internal_derive;
}

void
Lattice_element::set_needs_external_derive(bool value)
{
    needs_external_derive = value;
}

bool
Lattice_element::get_needs_external_derive() const
{
    return needs_external_derive;
}

double
Lattice_element::get_length() const
{
    std::map<std::string, double >::const_iterator iter =
            double_attributes.find(length_attribute_name);
    double retval = 0.0;
    if (iter != double_attributes.end()) {
        retval = iter->second;
    }
    return retval;
}

double
Lattice_element::get_bend_angle() const
{
    std::map<std::string, double >::const_iterator iter =
            double_attributes.find(bend_angle_attribute_name);
    double retval = 0.0;
    if (iter != double_attributes.end()) {
        retval = iter->second;
    }
    return retval;
}

long int
Lattice_element::get_revision() const
{
    return revision;
}

bool
Lattice_element::has_lattice() const
{
    return (lattice_ptr != 0);
}

void
Lattice_element::set_lattice(Lattice &lattice)
{
    lattice_ptr = &lattice;
}

Lattice &
Lattice_element::get_lattice()
{
    if(! has_lattice()) {
        throw std::runtime_error(
                    "Lattice_element::get_lattice: element not part of any lattice");
    }
    return *lattice_ptr;
}

Lattice const&
Lattice_element::get_lattice() const
{
    if(! has_lattice()) {
        throw std::runtime_error(
                    "Lattice_element::get_lattice: element not part of any lattice");
    }
    return *lattice_ptr;
}

