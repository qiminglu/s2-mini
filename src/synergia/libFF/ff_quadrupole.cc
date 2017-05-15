#include "ff_quadrupole.h"
#include "ff_algorithm.h"
#include "synergia/utils/gsvector.h"
#include "synergia/utils/logger.h"

double FF_quadrupole::get_reference_cdt(double length, double * k,
                                        Reference_particle &reference_particle) {
    double reference_cdt;
    if (length == 0) {
        reference_cdt = 0.0;
    } else {
        double reference_momentum = reference_particle.get_momentum();
        double m = reference_particle.get_mass();
        double step_length = length/steps;
        double step_strength[2] = { k[0]*step_length, k[1]*step_length };

        double x(reference_particle.get_state()[Bunch::x]);
        double xp(reference_particle.get_state()[Bunch::xp]);
        double y(reference_particle.get_state()[Bunch::y]);
        double yp(reference_particle.get_state()[Bunch::yp]);
        double cdt(reference_particle.get_state()[Bunch::cdt]);
        double dpop(reference_particle.get_state()[Bunch::dpop]);

        double cdt_orig = cdt;
        FF_algorithm::yoshida4<double, FF_algorithm::thin_quadrupole_unit<double>, 1 >
                ( x, xp, y, yp, cdt, dpop,
                  reference_momentum, m,
                  0.0,
                  step_length, step_strength, steps );
        reference_cdt = cdt - cdt_orig;
    }
    return reference_cdt;
}

void FF_quadrupole::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    throw std::runtime_error("Propagate JetParticle through a quadrupole element is yet to be implemented");

#if 0
    double length = slice.get_right() - slice.get_left();

    double k[2];
    k[0] = slice.get_lattice_element().get_double_attribute("k1", 0.0);
    k[1] = slice.get_lattice_element().get_double_attribute("k1s", 0.0);

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
    double reference_cdt = get_reference_cdt(length, k, reference_particle);
    double step_reference_cdt = reference_cdt/steps;
    double step_length = length/steps;
    double step_strength[2] = { k[0]*step_length, k[1]*step_length };
    double kl[2] = { k[0]*length, k[1]*length };

    if (length == 0.0) {
        FF_algorithm::thin_quadrupole_unit(x, xp, y, yp, kl);
    } else {
        FF_algorithm::yoshida4<TJet<double>, FF_algorithm::thin_quadrupole_unit<TJet<double> >, 1 >
                ( x, xp, y, yp, cdt, dpop,
                  reference_momentum, m,
                  step_reference_cdt,
                  step_length, step_strength, steps );
    }
    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, reference_momentum, m,
               reference_cdt);
#endif
}

void FF_quadrupole::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    // length
    double length = slice.get_right() - slice.get_left();

    // strength
    double k[2];
    k[0] = slice.get_lattice_element().get_double_attribute("k1", 0.0);
    k[1] = slice.get_lattice_element().get_double_attribute("k1s", 0.0);

    // offsets
    const double xoff = slice.get_lattice_element().get_double_attribute("hoffset", 0.0);
    const double yoff = slice.get_lattice_element().get_double_attribute("voffset", 0.0);

    const GSVector vxoff(xoff);
    const GSVector vyoff(yoff);

    // tilting
    double tilt = slice.get_lattice_element().get_double_attribute("tilt", 0.0);
    if (tilt != 0.0)
    {
        std::complex<double> ck2(k[0], -k[1]);
        ck2 = ck2 * exp(std::complex<double>(0.0, -2.0*tilt));
        k[0] = ck2.real();
        k[1] = ck2.imag();
    }

    // scaling
    Reference_particle const & ref_l = get_ref_particle_from_slice(slice);
    Reference_particle const & ref_b = bunch.get_reference_particle();

    double brho_l = ref_l.get_momentum() / ref_l.get_charge();  // GV/c
    double brho_b = ref_b.get_momentum() * (1.0 + ref_b.get_state()[Bunch::dpop]) / ref_l.get_charge();  // GV/c

    double scale = brho_l / brho_b;

    k[0] *= scale;
    k[1] *= scale;


    int local_num = bunch.get_local_num();
    MArray2d_ref particles = bunch.get_local_particles();

    double * xa, * xpa;
    double * ya, * ypa;
    double * cdta, * dpopa;

    bunch.set_arrays(xa, xpa, ya, ypa, cdta, dpopa);

    const int gsvsize = GSVector::size();
    const int num_blocks = local_num / gsvsize;
    const int block_last = num_blocks * gsvsize;

    if (length == 0.0) 
    {
#if 0
        #pragma omp parallel for
        #pragma acc kernels
        for (int part = 0; part < block_last; part += gsvsize) 
        {
            GSVector  x( &xa[part]);
            GSVector xp(&xpa[part]);
            GSVector  y( &ya[part]);
            GSVector yp(&ypa[part]);

            x -= vxoff;
            y -= vyoff;

            FF_algorithm::thin_quadrupole_unit(x, xp, y, yp, k);

            xp.store(&xpa[part]);
            yp.store(&ypa[part]);
        }
#endif

        //for (int part = block_last; part < local_num; ++part) 
        #pragma acc kernels
        for (int part = 0; part < local_num; ++part) 
        {
            double  x( xa[part]);
            double xp(xpa[part]);
            double  y( ya[part]);
            double yp(ypa[part]);

            x -= xoff;
            y -= yoff;

            FF_algorithm::thin_quadrupole_unit(x, xp, y, yp, k);

            xpa[part] = xp;
            ypa[part] = yp;
        }
    } 
    else 
    {
        // yoshida steps
        steps = (int)slice.get_lattice_element().get_double_attribute("yoshida_steps", 4.0);

        // params
        double reference_momentum = bunch.get_reference_particle().get_momentum();
        double m = bunch.get_mass();
        double reference_cdt = get_reference_cdt(length, k, bunch.get_reference_particle());
        double step_reference_cdt = reference_cdt/steps;
        double step_length = length/steps;
        double step_strength[2] = { k[0]*step_length, k[1]*step_length };

#if 0
        #pragma omp parallel for
        #pragma acc kernels
        for (int part = 0; part < block_last; part += gsvsize) 
        {
            GSVector    x(   &xa[part]);
            GSVector   xp(  &xpa[part]);
            GSVector    y(   &ya[part]);
            GSVector   yp(  &ypa[part]);
            GSVector  cdt( &cdta[part]);
            GSVector dpop(&dpopa[part]);

            x -= vxoff;
            y -= vyoff;

#if 0
            FF_algorithm::yoshida<GSVector, FF_algorithm::thin_quadrupole_unit<GSVector>, 6/*order*/, 1/*components*/ >
                    ( x, xp, y, yp, cdt, dpop,
                      reference_momentum, m,
                      step_reference_cdt,
                      step_length, step_strength, steps );
#endif

#if 1
            FF_algorithm::yoshida6<GSVector, FF_algorithm::thin_quadrupole_unit<GSVector>, 1 >
                    ( x, xp, y, yp, cdt, dpop,
                      reference_momentum, m,
                      step_reference_cdt,
                      step_length, step_strength, steps );
#endif

#if 0
            FF_algorithm::yoshida4<GSVector, FF_algorithm::thin_quadrupole_unit<GSVector>, 1 >
                    ( x, xp, y, yp, cdt, dpop,
                      reference_momentum, m,
                      step_reference_cdt,
                      step_length, step_strength, steps );
#endif

#if 0
            FF_algorithm::quadrupole_chef<GSVector>
                    ( x, xp, y, yp, cdt, dpop,
                      reference_momentum, m,
                      reference_cdt, length, k, 40 /* kicks */ );
#endif

            x += vxoff;
            y += vyoff;

               x.store(&xa[part]);
              xp.store(&xpa[part]);
               y.store(&ya[part]);
              yp.store(&ypa[part]);
             cdt.store(&cdta[part]);
            dpop.store(&dpopa[part]);
        }
#endif


        #pragma acc data copy(xa[0:local_num], xpa[0:local_num], ya[0:local_num], ypa[0:local_num], cdta[0:local_num], dpopa[0:local_num])
        {

        //for (int part = block_last; part < local_num; ++part) 
        #pragma acc kernels loop
        for (int part = 0; part < local_num; ++part) 
        {
            double    x(   xa[part]);
            double   xp(  xpa[part]);
            double    y(   ya[part]);
            double   yp(  ypa[part]);
            double  cdt( cdta[part]);
            double dpop(dpopa[part]);

            x -= xoff;
            y -= yoff;

            FF_algorithm::yoshida6<double, FF_algorithm::thin_quadrupole_unit<double>, 1 >
                    ( x, xp, y, yp, cdt, dpop,
                      reference_momentum, m,
                      step_reference_cdt,
                      step_length, step_strength, steps );

            x += xoff;
            y += yoff;

               xa[part] = x;
              xpa[part] = xp;
               ya[part] = y;
              ypa[part] = yp;
             cdta[part] = cdt;
            dpopa[part] = dpop;
        }
        }

        bunch.get_reference_particle().increment_trajectory(length);
    }
}

FF_quadrupole::~FF_quadrupole()
{

}

