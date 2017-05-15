#include "ff_drift.h"
#include "synergia/utils/gsvector.h"
#include "synergia/utils/logger.h"

FF_drift::FF_drift()
{
}

double FF_drift::get_reference_cdt(double length, Reference_particle & reference_particle)
{
    double x(reference_particle.get_state()[Bunch::x]);
    double xp(reference_particle.get_state()[Bunch::xp]);
    double y(reference_particle.get_state()[Bunch::y]);
    double yp(reference_particle.get_state()[Bunch::yp]);
    double cdt(reference_particle.get_state()[Bunch::cdt]);
    double dpop(reference_particle.get_state()[Bunch::dpop]);
    double reference_momentum = reference_particle.get_momentum();
    double m = reference_particle.get_mass();

    double cdt_orig = cdt;
    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, reference_momentum, m, 0.0);

    return cdt - cdt_orig;
}

void FF_drift::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    throw std::runtime_error("Propagate JetParticle through a drift element is yet to be implemented");

#if 0
    double length = slice.get_right() - slice.get_left();

    typedef PropagatorTraits<JetParticle>::State_t State_t;
    typedef PropagatorTraits<JetParticle>::Component_t Component_t;

    State_t& state = jet_particle.State();

    Component_t & x(state[Chef::x]);
    Component_t & xp(state[Chef::xp]);
    Component_t & y(state[Chef::y]);
    Component_t & yp(state[Chef::yp]);
    Component_t & cdt(state[Chef::cdt]);
    Component_t & dpop(state[Chef::dpop]);

    double reference_momentum = jet_particle.ReferenceMomentum();
    double m = jet_particle.Mass();

    Particle chef_particle(jet_particle);
    Reference_particle reference_particle(
                chef_particle_to_reference_particle(chef_particle));
    double reference_cdt = get_reference_cdt(length, reference_particle);

    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, reference_momentum, m,
               reference_cdt);
#endif
}

void FF_drift::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    const double  length = slice.get_right() - slice.get_left();
    const int  local_num = bunch.get_local_num();
    const double    mass = bunch.get_mass();
    const double ref_cdt = get_reference_cdt(length, bunch.get_reference_particle());

    Reference_particle & ref_b = bunch.get_reference_particle();
    const double   ref_p = ref_b.get_momentum() * (1.0 + ref_b.get_state()[Bunch::dpop]);

    double * xa, * xpa;
    double * ya, * ypa;
    double * cdta, * dpopa;

    bunch.set_arrays(xa, xpa, ya, ypa, cdta, dpopa);

    const int gsvsize = GSVector::size();
    const int num_blocks = local_num / gsvsize;
    const int block_last = num_blocks * gsvsize;

#if 0
    #pragma omp parallel for
    for (int part = 0; part < block_last; part += gsvsize) 
    {
        GSVector x(&xa[part]);
        GSVector xp(&xpa[part]);
        GSVector y(&ya[part]);
        GSVector yp(&ypa[part]);
        GSVector cdt(&cdta[part]);
        GSVector dpop(&dpopa[part]);

        FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, ref_p, mass, ref_cdt);

        x.store(&xa[part]);
        y.store(&ya[part]);
        cdt.store(&cdta[part]);
    }
#endif

    //for (int part = block_last; part < local_num; ++part) 
    #pragma omp parallel for
    #pragma acc data copy(xa[0:local_num], xpa[0:local_num], ya[0:local_num], ypa[0:local_num], cdta[0:local_num], dpopa[0:local_num])
    #pragma acc kernels
    for (int part = 0; part < local_num; ++part) 
    {
#if 0
        double x(xa[part]);
        double xp(xpa[part]);
        double y(ya[part]);
        double yp(ypa[part]);
        double cdt(cdta[part]);
        double dpop(dpopa[part]);

        FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, ref_p, mass, ref_cdt);

        xa[part] = x;
        ya[part] = y;
        cdta[part] = cdt;
#endif

        FF_algorithm::drift_unit(xa[part], xpa[part], ya[part], ypa[part], cdta[part], dpopa[part], length, ref_p, mass, ref_cdt);
    }

    bunch.get_reference_particle().increment_trajectory(length);
}

FF_drift::~FF_drift()
{

}
