#include "reference_particle.h"
#include "synergia/utils/floating_point.h"

Reference_particle::Reference_particle() :
    charge(0), four_momentum(1.0, 1.0), state(6),
            repetition(0), s(0), s_n(0)
{

}

Reference_particle::Reference_particle(int charge, double mass,
        double total_energy) :
    charge(charge), four_momentum(mass, total_energy),
            state(6, 0), repetition(0), s(0), s_n(0)
{
}

Reference_particle::Reference_particle(int charge,
        Four_momentum const & four_momentum_in) :
    charge(charge), four_momentum(four_momentum_in), state(6, 0),
            repetition(0), s(0), s_n(0)
{
}

Reference_particle::Reference_particle(int charge,
        Four_momentum const & four_momentum_in, Const_MArray1d_ref state) :
    charge(charge), four_momentum(four_momentum_in), state(state),
            repetition(0), s(0), s_n(0)
{
}

void
Reference_particle::set_four_momentum(Four_momentum const & four_momentum)
{
    this->four_momentum = four_momentum;
}

void
Reference_particle::set_state(Const_MArray1d_ref state)
{
    this->state = state;
}

void
Reference_particle::set_state(double x, double xp, double y, double yp, double cdt, double dpop)
{
    this->state[0] = x;
    this->state[1] = xp;
    this->state[2] = y;
    this->state[3] = yp;
    this->state[4] = cdt;
    this->state[5] = dpop;
}

void
Reference_particle::set_total_energy(double total_energy)
{
    four_momentum.set_total_energy(total_energy);
}

void
Reference_particle::increment_trajectory(double length)
{
    s_n += length;
}

void
Reference_particle::start_repetition()
{
    if (s == 0.0) {
        s = s_n;
    }
    if (s_n > 0.0) {
        repetition += 1;
    }
    s_n = 0.0;
}

void
Reference_particle::set_trajectory(int repetition, double repetition_length,
        double s)
{
    this->repetition = repetition;
    this->s = repetition_length;
    this->s_n = s;
}

int
Reference_particle::get_charge() const
{
    return charge;
}

double
Reference_particle::get_mass() const
{
    return four_momentum.get_mass();
}

Four_momentum const &
Reference_particle::get_four_momentum() const
{
    return four_momentum;
}

Const_MArray1d_ref
Reference_particle::get_state() const
{
    return state;
}

double
Reference_particle::get_gamma() const
{
    return four_momentum.get_gamma();
}

double
Reference_particle::get_beta() const
{
    return four_momentum.get_beta();
}

double
Reference_particle::get_momentum() const
{
    return four_momentum.get_momentum();
}

double
Reference_particle::get_total_energy() const
{
    return four_momentum.get_total_energy();
}

double
Reference_particle::get_s() const
{
    return repetition * s + s_n;
}

double
Reference_particle::get_s_n() const
{
    return s_n;
}

int
Reference_particle::get_repetition() const
{
    return repetition;
}

double
Reference_particle::get_repetition_length() const
{
    return s;
}

bool
Reference_particle::equal(Reference_particle const& reference_particle,
        double tolerance) const
{
    if (charge != reference_particle.get_charge()) {
        return false;
    }
    if (!four_momentum.equal(reference_particle.get_four_momentum(), tolerance)) {
        return false;
    }
    for (int i = 0; i < 6; ++i) {
        if (!floating_point_equal(state[i], reference_particle.get_state()[i],
                tolerance)) {
            return false;
        }
    }
    return true;
}
