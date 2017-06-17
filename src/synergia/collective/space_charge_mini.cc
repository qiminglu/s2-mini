
#include "space_charge_mini.h"

#include <cstring>
#include <omp.h>


inline int  fast_int_floor(double x)
{
    int ix = static_cast<int>(x);
    return x > 0.0 ? ix : ((x - ix == 0) ? ix : ix - 1);
}

inline bool ingrid(int x, int y, int z, int gx, int gy, int gz)
{
  return x>=0 && y>=0 && z>=0 && x<gx && y<gy && z<gz;
}

Space_charge_mini::Space_charge_mini(std::vector<int> const & grid_shape) 
: Collective_operator("space_charge_mini") 
, shape(grid_shape)
, rho(shape[0] * shape[1] * shape[2])
{

}

void Space_charge_mini::apply(Bunch & bunch, double time_step, Step & step, int verbosity, Logger & logger)
{
    get_charge_density_reduce(bunch, logger);
    logger << rho[0] << ", " << rho[1] << "\n";
}

void Space_charge_mini::get_charge_density(Bunch & bunch, Logger & logger)
{
    const int npart = bunch.get_local_num();
    svec<double> const & parts = bunch.get_local_particles();

    double * prho = rho.origin();

    std::vector<double> offs = { 0, 0, 0 };
    std::vector<double> size = { 2, 2, 2 };
    std::vector<double> h = { size[0]/shape[0], size[1]/shape[1], size[2]/shape[2] };

    double w0 = 1e4 * bunch.get_particle_charge() * 1.602176565e-19 / (h[0]*h[1]*h[2]);

    int gx = shape[0];
    int gy = shape[1];
    int gz = shape[2];

    int nc = gx * gy * gz;

    static int nt = 0;

    static double * rl = 0;
    static double * rs = 0;

    if (nt == 0)
    {
        #pragma omp parallel
        { nt = omp_get_num_threads(); }

        rl = new double[nt * nc];
        rs = new double[nc];
    }

    //logger << "nt = " << nt << "\n";

    double lx = offs[0] - size[0] / 2.0;
    double ly = offs[1] - size[1] / 2.0;
    double lz = offs[2] - size[2] / 2.0;

    double cx = h[0];
    double cy = h[1];
    double cz = h[2];

#if 0
    #pragma omp parallel for reduction(+:S[:10])
    for (int i=0; i<npart; ++i)
    {
        int ix, iy, iz;
        double ox, oy, oz;

        double x = parts[0 * npart +i];
        double y = parts[2 * npart +i];
        double z = parts[4 * npart +i];

        scaled_location = (x - lx) / cx - 0.5;
        ix = fast_int_floor(scaled_location);
        ox = scaled_location - ix;

        scaled_location = (y - ly) / cy - 0.5;
        iy = fast_int_floor(scaled_location);
        oy = scaled_location - iy;

        scaled_location = (z - lz) / cz - 0.5;
        iz = fast_int_floor(scaled_location);
        oz = scaled_location - iz;
    }
#endif

    #pragma omp parallel shared(parts, lx, ly, lz, cx, cy, cz, w0, gx, gy, gz, rl, rs, nc, prho)
    {
        int nt = omp_get_num_threads();
        int it = omp_get_thread_num();

        int np = npart;
        int le = npart / nt;
        int ps = it * le;
        int pe = (it == nt-1) ? np : (it+1)*le;

        // zero the worksheet
        std::memset(rl + it*nc, 0, sizeof(double)*nc);

        double scaled_location;
        double ox, oy, oz, w;
        int    ix, iy, iz;

        for(int n=ps; n<pe; ++n)
        {
            scaled_location = (parts[0*npart + n] - lx) / cx - 0.5;
            ix = fast_int_floor(scaled_location);
            ox = scaled_location - ix;

            scaled_location = (parts[2*npart + n] - ly) / cy - 0.5;
            iy = fast_int_floor(scaled_location);
            oy = scaled_location - iy;

            scaled_location = (parts[4*npart + n] - lz) / cz - 0.5;
            iz = fast_int_floor(scaled_location);
            oz = scaled_location - iz;

            int base = iz * gx * gy;

#if 0
            if( ingrid(ix  , iy  , iz, gx, gy, gz) ) rl[(base + (iy  )*gx + ix  )*nt+it] += w0*(1-ox)*(1-oy)*(1-oz);
            if( ingrid(ix+1, iy  , iz, gx, gy, gz) ) rl[(base + (iy  )*gx + ix+1)*nt+it] += w0*(  ox)*(1-oy)*(1-oz);
            if( ingrid(ix  , iy+1, iz, gx, gy, gz) ) rl[(base + (iy+1)*gx + ix  )*nt+it] += w0*(1-ox)*(  oy)*(1-oz);
            if( ingrid(ix+1, iy+1, iz, gx, gy, gz) ) rl[(base + (iy+1)*gx + ix+1)*nt+it] += w0*(  ox)*(  oy)*(1-oz);
#endif
            if( ingrid(ix  , iy  , iz, gx, gy, gz) ) rl[it*nt + (base + (iy  )*gx + ix  )] += w0*(1-ox)*(1-oy)*(1-oz);
            if( ingrid(ix+1, iy  , iz, gx, gy, gz) ) rl[it*nt + (base + (iy  )*gx + ix+1)] += w0*(  ox)*(1-oy)*(1-oz);
            if( ingrid(ix  , iy+1, iz, gx, gy, gz) ) rl[it*nt + (base + (iy+1)*gx + ix  )] += w0*(1-ox)*(  oy)*(1-oz);
            if( ingrid(ix+1, iy+1, iz, gx, gy, gz) ) rl[it*nt + (base + (iy+1)*gx + ix+1)] += w0*(  ox)*(  oy)*(1-oz);

            base = (iz+1) * gx * gy;

#if 0
            if( ingrid(ix  , iy  , iz+1, gx, gy, gz) ) rl[(base + (iy  )*gx + ix  )*nt+it] += w0*(1-ox)*(1-oy)*(oz);
            if( ingrid(ix+1, iy  , iz+1, gx, gy, gz) ) rl[(base + (iy  )*gx + ix+1)*nt+it] += w0*(  ox)*(1-oy)*(oz);
            if( ingrid(ix  , iy+1, iz+1, gx, gy, gz) ) rl[(base + (iy+1)*gx + ix  )*nt+it] += w0*(1-ox)*(  oy)*(oz); 
            if( ingrid(ix+1, iy+1, iz+1, gx, gy, gz) ) rl[(base + (iy+1)*gx + ix+1)*nt+it] += w0*(  ox)*(  oy)*(oz); 
#endif
            if( ingrid(ix  , iy  , iz+1, gx, gy, gz) ) rl[it*nt + (base + (iy  )*gx + ix  )] += w0*(1-ox)*(1-oy)*(oz);
            if( ingrid(ix+1, iy  , iz+1, gx, gy, gz) ) rl[it*nt + (base + (iy  )*gx + ix+1)] += w0*(  ox)*(1-oy)*(oz);
            if( ingrid(ix  , iy+1, iz+1, gx, gy, gz) ) rl[it*nt + (base + (iy+1)*gx + ix  )] += w0*(1-ox)*(  oy)*(oz); 
            if( ingrid(ix+1, iy+1, iz+1, gx, gy, gz) ) rl[it*nt + (base + (iy+1)*gx + ix+1)] += w0*(  ox)*(  oy)*(oz); 
        }

        #pragma omp barrier

        // reduction
        le = gz/nt;
        ps = it*le;
        pe = (it == nt-1) ? gz : (ps + le);

        for(int z=ps; z<pe; ++z)
        {
            for(int y=0; y<gy; ++y)
            {
                for(int x=0; x<gx; ++x)
                {
                    w = 0.0;
                    for(int n=0; n<nt; ++n) w += rl[(z*gx*gy + y*gx + x)*nt+n];
                    prho[z*gx*gy + y*gx +x] = w;
                }
            }
        }

        #pragma omp barrier

    } //  end of #pragma parallel

    //memcpy( prho, rs, sizeof(double)*nc );
}
void Space_charge_mini::get_charge_density_reduce(Bunch & bunch, Logger & logger)
{
    const int npart = bunch.get_local_num();
    svec<double> const & parts = bunch.get_local_particles();

    double * prho = rho.origin();

    std::vector<double> offs = { 0, 0, 0 };
    std::vector<double> size = { 2, 2, 2 };
    std::vector<double> h = { size[0]/shape[0], size[1]/shape[1], size[2]/shape[2] };

    double w0 = 1e4 * bunch.get_particle_charge() * 1.602176565e-19 / (h[0]*h[1]*h[2]);

    int gx = shape[0];
    int gy = shape[1];
    int gz = shape[2];

    int nc = gx * gy * gz;

    //logger << "nt = " << nt << "\n";

    double lx = offs[0] - size[0] / 2.0;
    double ly = offs[1] - size[1] / 2.0;
    double lz = offs[2] - size[2] / 2.0;

    double cx = h[0];
    double cy = h[1];
    double cz = h[2];

    //#pragma omp parallel for reduction(+:prho[:10])
    #pragma omp parallel for shared(parts, gx, gy, gz, lx, ly, lz, cx, cy, cz) reduction(+:prho[:nc])
    for (int i=0; i<npart; ++i)
    {
        int ix, iy, iz;
        double ox, oy, oz;

        double x = parts[0 * npart +i];
        double y = parts[2 * npart +i];
        double z = parts[4 * npart +i];

        double scaled_location;

        scaled_location = (x - lx) / cx - 0.5;
        ix = fast_int_floor(scaled_location);
        ox = scaled_location - ix;

        scaled_location = (y - ly) / cy - 0.5;
        iy = fast_int_floor(scaled_location);
        oy = scaled_location - iy;

        scaled_location = (z - lz) / cz - 0.5;
        iz = fast_int_floor(scaled_location);
        oz = scaled_location - iz;

        int base = iz * gx * gy;

        if( ingrid(ix  , iy  , iz, gx, gy, gz) ) prho[base + (iy  )*gx + ix  ] += w0*(1-ox)*(1-oy)*(1-oz);
        if( ingrid(ix+1, iy  , iz, gx, gy, gz) ) prho[base + (iy  )*gx + ix+1] += w0*(  ox)*(1-oy)*(1-oz);
        if( ingrid(ix  , iy+1, iz, gx, gy, gz) ) prho[base + (iy+1)*gx + ix  ] += w0*(1-ox)*(  oy)*(1-oz);
        if( ingrid(ix+1, iy+1, iz, gx, gy, gz) ) prho[base + (iy+1)*gx + ix+1] += w0*(  ox)*(  oy)*(1-oz);

        base = (iz+1) * gx * gy;

        if( ingrid(ix  , iy  , iz+1, gx, gy, gz) ) prho[base + (iy  )*gx + ix  ] += w0*(1-ox)*(1-oy)*(oz);
        if( ingrid(ix+1, iy  , iz+1, gx, gy, gz) ) prho[base + (iy  )*gx + ix+1] += w0*(  ox)*(1-oy)*(oz);
        if( ingrid(ix  , iy+1, iz+1, gx, gy, gz) ) prho[base + (iy+1)*gx + ix  ] += w0*(1-ox)*(  oy)*(oz); 
        if( ingrid(ix+1, iy+1, iz+1, gx, gy, gz) ) prho[base + (iy+1)*gx + ix+1] += w0*(  ox)*(  oy)*(oz); 

    }
}
