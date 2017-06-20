
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
    //get_charge_density_reduce(bunch, logger);
    get_charge_density(bunch, logger);
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



void Space_charge_mini::get_charge_density_interleaved(Bunch & bunch, Logger & logger)
{
    const int npart = bunch.get_local_num();
    svec<double> const & parts = bunch.get_local_particles();

    double * prho = rho.origin();

    std::array<double, 3> offs = { 0, 0, 0 };
    std::array<double, 3> size = { 2, 2, 2 };
    std::array<double, 3> h = { size[0]/shape[0], size[1]/shape[1], size[2]/shape[2] };

    double w0 = 1e4 * bunch.get_particle_charge() * 1.602176565e-19 / (h[0]*h[1]*h[2]);

    int gx = shape[0];
    int gy = shape[1];
    int gz = shape[2];

    double * po = new double[npart*3];  // offx, offy, offz
    int    * pi = new int[npart];       // cell index of particles
    int    * pc = new int[gx*gy*gz+1];  // accumulated particle count in cells
    int    * count = new int[gx*gy*gz]; // particle count in cells
    int    * pll = new int[npart];      // indexed list

    static int nt = 0;

    if( nt==0 )
    {
        #pragma omp parallel
        { nt = omp_get_num_threads(); }
    }

    #pragma omp parallel for
    for(int i=0; i<gx*gy*gz+1; ++i) pc[i] = 0;
    //memset( pc   , 0, sizeof(int)*(gx*gy*gz+1) );

    #pragma omp parallel for
    for(int i=0; i<gx*gy*gz; ++i) count[i] = 0;
    //memset( count, 0, sizeof(int)*(gx*gy*gz) );

    // zero first
    #pragma omp parallel for
    for(int i=0; i<gx*gy*gz; ++i) prho[i] = 0.0;

    int ix, iy, iz;
    double offx, offy, offz;

    double lx = offs[2] - size[2] / 2.0;
    double ly = offs[1] - size[1] / 2.0;
    double lz = offs[0] - size[0] / 2.0;

    double cx = h[2];
    double cy = h[1];
    double cz = h[0];

   {
        // decl. for private variables
        int n, i, c_idx, x, y, z;
        int tid, nthreads, seg;
        double scaled_location;
        double w, ox, oy, oz;
      
        int * p;
        double * ws0;
        double * ws1;

        // count
        #pragma omp parallel for \
            shared( pi, po ) \
            private(n, c_idx, p, ix, iy, iz, offx, offy, offz, scaled_location)
        for(n=0; n<npart; ++n )
        {
            scaled_location = (parts[0*npart + n] - lx) / cx - 0.5;
            ix = fast_int_floor(scaled_location);
            offx = scaled_location - ix;

            scaled_location = (parts[2*npart + n] - ly) / cy - 0.5;
            iy = fast_int_floor(scaled_location);
            offy = scaled_location - iy;

            scaled_location = (parts[4*npart + n] - lz) / cz - 0.5;
            iz = fast_int_floor(scaled_location);
            offz = scaled_location - iz;

            if( ix<0 || iy<0 || iz<0 || ix>=gx || iy>=gy || iz>=gz )
            {
              pi[n] = -1;
            }
            else
            {
              c_idx = iz*gx*gy + iy*gx + ix;
              pi[n] = c_idx;

              po[n*3+0] = offx;
              po[n*3+1] = offy;
              po[n*3+2] = offz;

              // gcc build-in atomic add
              __sync_fetch_and_add(pc+c_idx+1, 1);
              
              // another choice is to use the omp atomic
              //#pragma omp atomic
              //++pc[c_idx+1];
            }
        }

        // accumulate particle counts
        for(i=1; i<=gx*gy*gz; ++i)
          pc[i] += pc[i-1];

        int idx, pos;

        // build indexed list in parallel
        #pragma omp parallel for shared(pi, count) private(i, idx, pos)
        for(i=0; i<npart; ++i)
        {
          idx = pi[i];
          if( idx==-1 ) continue;
          pos = pc[idx] + __sync_fetch_and_add( count+idx, 1 );
          pll[pos] = i;
        }

        // deposit to cells
        #pragma omp parallel \
            shared( rho, pi, po, pc, count, pll  ) \
            private(x, y, z, ox, oy, oz, n, i, c_idx, w, ws0, ws1, tid, nthreads, seg) 
        {
          tid = omp_get_thread_num();
          nthreads = omp_get_num_threads();

          // temp work table to hold intermediate results
          ws0 = new double[(gx+1)*(gy+1)]();
          ws1 = new double[(gx+1)*(gy+1)]();

          // interleaved partitioning along z-axis to have 
          // evenly balanced loads for each thread
          for(z=tid; z<gz; z+=nthreads)
          {
            memset(ws0, 0, sizeof(double)*(gx+1)*(gy+1));
            memset(ws1, 0, sizeof(double)*(gx+1)*(gy+1));

            for(x=0; x<gx; ++x)
            {
              for(y=0; y<gy; ++y)
              {
                c_idx = z*gx*gy + y*gx + x;

                for(n=0; n<count[c_idx]; ++n)
                {
                  i = pll[pc[c_idx]+n]; // index of the particle

                  ox = po[i*3+0]; oy = po[i*3+1]; oz = po[i*3+2];

                  w = w0 * (1-ox) * (1-oy) * (1-oz); ws0[(y  )*(gx+1) + x  ] += w;
                  w = w0 * (  ox) * (1-oy) * (1-oz); ws0[(y  )*(gx+1) + x+1] += w;
                  w = w0 * (1-ox) * (  oy) * (1-oz); ws0[(y+1)*(gx+1) + x  ] += w;
                  w = w0 * (  ox) * (  oy) * (1-oz); ws0[(y+1)*(gx+1) + x+1] += w;
                
                  w = w0 * (1-ox) * (1-oy) * (  oz); ws1[(y  )*(gx+1) + x  ] += w;
                  w = w0 * (  ox) * (1-oy) * (  oz); ws1[(y  )*(gx+1) + x+1] += w;
                  w = w0 * (1-ox) * (  oy) * (  oz); ws1[(y+1)*(gx+1) + x  ] += w;
                  w = w0 * (  ox) * (  oy) * (  oz); ws1[(y+1)*(gx+1) + x+1] += w;
                } 

              } //end of y
            } // end of x
        
            // write to global memory
            for(x=0; x<gx; ++x)
              for(y=0; y<gy; ++y)
                prho[z*gx*gy + y*gx +x] += ws0[ y*(gx+1) + x];
                //rho[z][y][x] += ws0[ y*(gx+1) + x ];

            // a synchronization is needed to avoid data pollution
            #pragma omp barrier

            if( z!=gz-1 )
              for(x=0; x<gx; ++x)
                for(y=0; y<gy; ++y)
                  prho[(z+1)*gx*gy + y*gx +x] += ws1[ y*(gx+1) + x];
                  //rho[z+1][y][x] += ws1[ y*(gx+1) + x ];

          } // end of z loop


          delete [] ws0;
          delete [] ws1;
        } 

    }

    delete [] po;
    delete [] pi;
    delete [] pc;
    delete [] count;
    delete [] pll;
}



