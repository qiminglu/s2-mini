#include "ff_drift.h"
#include "synergia/utils/gsvector.h"
#include "synergia/utils/logger.h"


inline void drift_unit(double *x, double xp, double *y, double yp, double *cdt, double dpop, double length, double ref_p, double mass, double ref_cdt)
{
        double dp = dpop + 1.0;
        double inv_npz = 1.0 / /*sqrt*/(dp * dp - xp * xp - yp * yp);
        double lxpr = xp * length * inv_npz;
        double lypr = yp * length * inv_npz;
        double D2 = lxpr * lxpr + length * length + lypr * lypr;
        double p = dp * ref_p;
        double E2 = 1.0;//p * p + mass * mass;
        double ibeta2 = E2 / (p * p);
        *x = *x + lxpr;
        *y = *y + lypr;
        //cdt += sqrt(D2 / beta2) - reference_cdt;
        //cdt = cdt + sig * sqrt(D2 * ibeta2) - vrc;
        //*cdt = *cdt + D2 * ibeta2 - ref_cdt;
}

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
          double  length = slice.get_right() - slice.get_left();
    const int  local_num = bunch.get_local_num();
    const double    mass = bunch.get_mass();
    const double ref_cdt = get_reference_cdt(length, bunch.get_reference_particle());

    Reference_particle & ref_b = bunch.get_reference_particle();
    const double   ref_p = ref_b.get_momentum() * (1.0 + ref_b.get_state()[Bunch::dpop]);

    double * xa, * xpa;
    double * ya, * ypa;
    double * cdta, * dpopa;

    bunch.set_arrays(xa, xpa, ya, ypa, cdta, dpopa);

    double * __restrict parts = bunch.get_local_particles().origin();

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
    //#pragma acc data copy(xa[0:local_num], ya[0:local_num], cdta[0:local_num]), copyin(xpa[0:local_num], ypa[0:local_num], dpopa[0:local_num], length, ref_p, mass, ref_cdt)
    #pragma acc parallel loop present(parts)
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

        FF_algorithm::drift_unit(
                parts[local_num*0 + part], 
                parts[local_num*1 + part], 
                parts[local_num*2 + part], 
                parts[local_num*3 + part], 
                parts[local_num*4 + part], 
                parts[local_num*5 + part], 
                length, ref_p, mass, ref_cdt );

        //FF_algorithm::drift_unit(xa[part], xpa[part], ya[part], ypa[part], cdta[part], dpopa[part], length, ref_p, mass, ref_cdt);
        //FF_algorithm::drift_unit_d(&xa[part], xpa[part], &ya[part], ypa[part], &cdta[part], dpopa[part], length, ref_p, mass, ref_cdt);

        //drift_unit(&xa[part], xpa[part], &ya[part], ypa[part], &cdta[part], dpopa[part], length, ref_p, mass, ref_cdt);
#if 0
        double dp = dpopa[part] + 1.0;
        double inv_npz = 1.0 / /*sqrt*/(dp * dp - xpa[part] * xpa[part] - ypa[part] * ypa[part]);
        double lxpr = xpa[part] * length * inv_npz;
        double lypr = ypa[part] * length * inv_npz;
        double D2 = lxpr * lxpr + length * length + lypr * lypr;
        double p = dp * ref_p;
        double E2 = p * p + mass * mass;
        //double beta2 = p*p / E2;
        double ibeta2 = E2 / (p * p);
        xa[part] = xa[part] + lxpr;
        ya[part] = ya[part] + lypr;
        //cdt += sqrt(D2 / beta2) - reference_cdt;
        //cdt = cdt + sig * sqrt(D2 * ibeta2) - vrc;
        cdta[part] = cdta[part] + D2 * ibeta2 - ref_cdt;
#endif
    }

    bunch.get_reference_particle().increment_trajectory(length);
}

FF_drift::~FF_drift()
{

}
