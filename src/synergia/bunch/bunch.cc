#include "bunch.h"
#include "synergia/utils/parallel_utils.h"

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <sstream>

const int Bunch::x;
const int Bunch::xp;
const int Bunch::y;
const int Bunch::yp;
const int Bunch::z;
const int Bunch::zp;
const int Bunch::cdt;
const int Bunch::dpop;
const int Bunch::id;


class Particle_id_offset
{
private:
    int offset;
public:
    Particle_id_offset() :
        offset(0)
    {
    }

    int
    get(int request_num, Commxx const & comm)
    {
        MPI_Bcast((void *) &offset, 1, MPI_INT, 0, comm.get());
        int old_offset = offset;
        int total_num;
        MPI_Reduce((void*) &request_num, (void*) &total_num, 1, MPI_INT,
                MPI_SUM, 0, comm.get());
        offset += total_num;
        return old_offset;
    }

};


static Particle_id_offset particle_id_offset;

void
Bunch::assign_ids(int local_offset)
{
    int global_offset, request_num;
    if (comm_sptr->get_rank() == 0) {
        request_num = total_num;
    } else {
        request_num = 0;
    }
    global_offset = particle_id_offset.get(request_num, *comm_sptr);
    for (int i = 0; i < local_num; ++i) {
        local_particles[local_num * id + i] = i + local_offset + global_offset;
    }
}



template<class T, size_t C, int I>
    struct Sortable2d
    {
        typedef T data_type;
    };

template<class T, size_t C, int I>
    struct Sortable2d<T*, C, I >
    {
        typedef T data_type;
        typedef T arr_type[][C];
        typedef T row_type[C];
        struct Row
        {
            row_type data;
        };
        typedef Row cols_type[];

        Sortable2d(double* t, size_t sz) :
            ptr_(t), rows_(sz)
        {
        }

        struct Less
        {
            bool
            operator()(Row const& a, Row const& b)
            {
                return a.data[I] < b.data[I];
            }
        };

        Row*
        begin()
        {
            return (Row*) ptr_;
        }

        Row*
        end()
        {
            return (Row*) (ptr_ + (rows_ * C));
        }

        double* ptr_;
        size_t rows_;
    };

void
Bunch::construct(int particle_charge, int total_num, double real_num)
{
    sort_counter = 0;
    sort_period = 10000;

    this->particle_charge = particle_charge;
    this->total_num = total_num;
    this->real_num = real_num;

    state = fixed_z;
    converter_ptr = &default_converter;

    if (comm_sptr->has_this_rank()) 
    {
        local_num = decompose_1d_local(*comm_sptr, total_num);
        std::vector<int > offsets(comm_sptr->get_size()), counts(comm_sptr->get_size());
        decompose_1d(*comm_sptr, total_num, offsets, counts);
        local_num = counts[comm_sptr->get_rank()];

        if (local_num != total_num)
        {
            local_particles.resize( local_num * 7 );
            std::cout << "resize local particles\n";
        }

        #pragma omp parallel for
        for (int i=0; i<local_num; ++i)
        {
            for(int j=0; j<7; ++j)
            {
#if 0
                (*local_particles)[i][j] = 0.0;
#endif
            }
        }

        assign_ids(offsets[comm_sptr->get_rank()]);
    } 
    else 
    {
        local_num = 0;
    }
}

Bunch::Bunch(Reference_particle const& reference_particle, int total_num, double real_num, Commxx_sptr comm_sptr) 
: z_period_length(0.0)
, z_periodic(0)
, reference_particle(reference_particle)
, local_particles(total_num * 7)
, local_num(total_num)
, real_num(real_num)
, bucket_index(0)
, comm_sptr(comm_sptr)
, default_converter()
{
    construct(reference_particle.get_charge(), total_num, real_num);
}

void
Bunch::set_particle_charge(int particle_charge)
{
    this->particle_charge = particle_charge;
}

void
Bunch::set_real_num(double real_num)
{
    this->real_num = real_num;
}

void
Bunch::set_local_num(int local_num)
{
#if 0
    if (local_num > this->local_num) 
    {
        // throw std::runtime_error("set num not supported");

        double * prev_storage = storage;
        double * prev_alt_storage = alt_storage;

        MArray2d_ref * prev_local_particles = local_particles;
        MArray2d_ref * prev_alt_local_particles = alt_local_particles;

        storage = (double*)boost::alignment::aligned_alloc(8 * sizeof(double), local_num * 7 * sizeof(double));
        alt_storage = (double*)boost::alignment::aligned_alloc(8 * sizeof(double), local_num * 7 * sizeof(double));

        local_particles = new MArray2d_ref(storage, boost::extents[local_num][7], boost::fortran_storage_order());
        alt_local_particles = new MArray2d_ref(alt_storage, boost::extents[local_num][7], boost::fortran_storage_order());

        for (int i=0; i<this->local_num; ++i) {
            for (int j=0; j<7; ++j) {
                (*local_particles)[i][j] = (*prev_local_particles)[i][j];
                (*alt_local_particles)[i][j] = (*prev_alt_local_particles)[i][j];
            }
        }

        delete [] prev_storage;
        delete [] prev_alt_storage;

        delete prev_local_particles;
        delete prev_alt_local_particles;


#if 0
        MArray2d *prev_local_particles = local_particles;
        int prev_local_num = this->local_num;
         local_particles = new MArray2d(boost::extents[local_num][7],
                 boost::fortran_storage_order());
         alt_local_particles = new MArray2d(boost::extents[local_num][7],
                 boost::fortran_storage_order());
         (*local_particles)[ boost::indices[range(0,prev_local_num)][range()] ] =
                (*prev_local_particles)[ boost::indices[range(0,prev_local_num)][range()] ];
        delete prev_local_particles;
#endif
    }

    this->local_num = local_num;
#endif
}

void
Bunch::update_total_num()
{
    int old_total_num = total_num;
    MPI_Allreduce(&local_num, &total_num, 1, MPI_INT, MPI_SUM, comm_sptr->get());
    if (old_total_num != 0.0) {
        real_num = (total_num * real_num) / old_total_num;
    } else {
        real_num = 0.0;
    }
}

void
Bunch::set_sort_period(int period)
{
    sort_period = period;
    sort_counter = period;
}

namespace {
double * semi_global_t;
size_t semi_global_start_pos;

inline bool do_compare(unsigned int const& a, unsigned int const& b)
{
    bool retval = semi_global_t[semi_global_start_pos+a] <
            semi_global_t[semi_global_start_pos+b];
    return retval;
}

void do_sort(double * t, size_t rows, size_t cols, size_t ord_col)
{
    semi_global_t = t;
    std::vector<unsigned int> index(cols);
    // c++ 11
    // unsigned int ind=0;
    //generate(index.begin(),index.end(), [&]() { return ind++; });
    for(int i=0; i<cols; ++i) {
        index[i] = i;
    }
    semi_global_start_pos = ord_col*cols;
    std::sort(index.begin(),index.end(), &do_compare);

    // swap all values in each row according to index order
    for(size_t r=0; r<rows; ++r) {
        double *start = t+(r*cols), *end = t+((r+1)*cols);
        std::vector<double> temp(start, end);
        for(size_t i=0; i<index.size(); ++i) {
            start[i] = temp[index[i]];
        }
    }
}
}

void
Bunch::sort(int index)
{
#if 0
    if ((index<0) || (index>6)) {
        throw std::runtime_error("Bunch::sort: invalid index");
    }
    do_sort(local_particles->origin(), 7, local_num, index);
    sort_counter = sort_period;
#endif
}

void
Bunch::periodic_sort(int index)
{
    if (sort_counter == 0) {
        sort(index);
    } else {
        --sort_counter;
    }
}

void
Bunch::set_converter(Fixed_t_z_converter &converter)
{
    this->converter_ptr = &converter;
}

void
Bunch::convert_to_state(State state)
{
    if (this->state != state) {
        if (this->state == fixed_z_lab) {
            if (state == fixed_t_lab) {
                converter_ptr->from_z_lab_to_t_lab(*this);
            }
            else if ( state == fixed_t_bunch) {
                converter_ptr->from_z_lab_to_t_bunch(*this);
            }
            // else if ( state == fixed_z_bunch) {
            //    converter_ptr->from_z_lab_to_z_bunch(*this);
            //}
            else {
                std::cout<<" state to convert to="<<state<<std::endl;
                std::cout<<" initial state ="<<this->state<<std::endl;
                throw std::runtime_error("Unknown state in Bunch::convert_to_state, case 1");
            }
        }
        else if (this->state == fixed_z_bunch) {
            throw std::runtime_error("state z_bunch not implemented yet in Bunch::convert_to_state");
        }
        else if (this->state == fixed_t_lab) {
            if (state == fixed_z_lab ) {
                converter_ptr->from_t_lab_to_z_lab(*this);
            }
            //else if (state == fixed_z_bunch) {
            //    converter_ptr->from_t_lab_to_z_bunch(*this);
            //}
            else if (state == fixed_t_bunch) {
                converter_ptr->from_t_lab_to_t_bunch(*this);
            }
            else {
                std::cout<<" state to convert to="<<state<<std::endl;
                std::cout<<" initial state ="<<this->state<<std::endl;
                throw std::runtime_error("Unknown state in Bunch::convert_to_state, case 2");
            }
        }
        else if (this->state == fixed_t_bunch) {
            if (state == fixed_z_lab ) {
                converter_ptr->from_t_bunch_to_z_lab(*this);
            }
            //else if (state == fixed_z_bunch ) {
            //    converter_ptr->from_t_bunch_to_z_bunch(*this);
            //}
            else if (state == fixed_t_lab ) {
                converter_ptr->from_t_bunch_to_t_lab(*this);
            }
            else {
                std::cout<<" state to convert to="<<state<<std::endl;
                std::cout<<" initial state ="<<this->state<<std::endl;
                throw std::runtime_error("Unknown state in Bunch::convert_to_state, case 3");
            }
        }
    this->state =state;
    }

}



Reference_particle &
Bunch::get_reference_particle()
{
    return reference_particle;
}

Reference_particle const&
Bunch::get_reference_particle() const
{
    return reference_particle;
}

MArray2d_ref
Bunch::get_local_particles()
{
    return local_particles;
}

Const_MArray2d_ref
Bunch::get_local_particles() const
{
    return local_particles;
}

int
Bunch::get_particle_charge() const
{
    return particle_charge;
}

double
Bunch::get_mass() const
{
    return reference_particle.get_mass();
}

double
Bunch::get_real_num() const
{
    return real_num;
}


double
 Bunch::get_z_period_length() const
{
    return z_period_length;
}

bool
 Bunch::is_z_periodic() const
{
    return z_periodic;
}

int
Bunch::get_local_num() const
{
    return local_num;
}

int
Bunch::get_total_num() const
{
    return total_num;
}

int
Bunch::get_sort_period() const
{
    return sort_period;
}



void
Bunch::set_bucket_index(int index)
{
    this->bucket_index=index;
}

int
Bunch::get_bucket_index() const
{
return bucket_index;
}



Bunch::State
Bunch::get_state() const
{
    return state;
}

Commxx const&
Bunch::get_comm() const
{
    return *comm_sptr;
}

Commxx_sptr
Bunch::get_comm_sptr() const
{
    return comm_sptr;
}

void
Bunch::inject(Bunch const& bunch)
{
#if 0
    const double weight_tolerance = 1.0e-10;
    const double particle_tolerance = 1.0e-14;

    // The charge and mass of the bunch particles must match
    if (particle_charge != bunch.get_particle_charge()) {
        throw std::runtime_error(
                "Bunch.inject: bunch particle charges do not match.");
    }
    if (std::abs(reference_particle.get_mass()/
                 bunch.get_reference_particle().get_mass() - 1.0) > particle_tolerance) {
        throw std::runtime_error(
                "Bunch:inject: bunch particle masses do not match.");
    }
    // can only check particle weight if total_num is nonzero
    if (total_num == 0) {
        // target bunch is empty.  Set the weights from the injected bunch
        real_num = bunch.get_real_num();
        total_num = bunch.get_total_num();
    } else if (std::abs(real_num/total_num - bunch.get_real_num()/bunch.get_total_num())
        > weight_tolerance) {
        throw std::runtime_error(
                "Bunch.inject: macroparticle weight of injected bunch does not match.");
    }
    int old_local_num = local_num;
    set_local_num(old_local_num + bunch.get_local_num());
    Const_MArray2d_ref injected_particles(bunch.get_local_particles());
    double target_momentum = reference_particle.get_momentum();
    double injected_momentum = bunch.get_reference_particle().get_momentum();
    MArray1d ref_state_diff(boost::extents[6]);
    MArray1d target_state(boost::extents[6]);
    MArray1d injected_state(boost::extents[6]);

    for (int i = 0; i < 6; ++i) {
        ref_state_diff[i] = bunch.get_reference_particle().get_state()[i]
                - reference_particle.get_state()[i];
    }

    for (int i = 0; i < 6; ++i) {
        target_state[i] = reference_particle.get_state()[i];
        injected_state[i] = bunch.get_reference_particle().get_state()[i];
    }

    for (int part = 0; part < bunch.get_local_num(); ++part) {
        // space-like coordinates
        for (int i = 0; i < 6; i += 2) {
            (*local_particles)[old_local_num + part][i]
                    = injected_particles[part][i] + ref_state_diff[i];
        }

        // npx and npy coordinates are scaled with p_ref which can be different
        // for different bunches
        for (int i = 1; i < 4; i += 2) {
            (*local_particles)[old_local_num + part][i] =
                    (injected_momentum/target_momentum) *
                    (injected_particles[part][i] - injected_state[i]) + target_state[i];
        }

        // ndp coordinate is delta-p scaled with pref
        (*local_particles)[old_local_num + part][5] =
                (injected_momentum/target_momentum) *
                (1.0 + injected_particles[part][5] - injected_state[5]) + target_state[5] - 1.0;

        (*local_particles)[old_local_num + part][Bunch::id]
                = injected_particles[part][Bunch::id];
    }
    update_total_num();
#endif
}

void Bunch::check_pz2_positive()
{
#if 0
    if (this->state == fixed_z_lab) {
        int local_num = get_local_num();
        MArray2d_ref particles = get_local_particles();
        for (int part = 0; part < local_num; ++part) {
            double  pzop2=(1.+particles[part][5])*(1.+particles[part][5])-
                particles[part][1]*particles[part][1]-particles[part][3]*particles[part][3];
            if (pzop2<0.)  {
                std::cout<<"pzop^2="<<pzop2<<std::endl;
                throw std::runtime_error( " check pz2:  pz square cannot be negative!");
            }

        }
    }
#endif
}

void Bunch::set_arrays(double * &xa, double * &xpa,
                       double * &ya, double * &ypa,
                       double * &cdta, double * &dpopa)
{
    double *origin = &local_particles[0];
    int stride = local_num;
    xa    = origin + stride*Bunch::x;
    xpa   = origin + stride*Bunch::xp;
    ya    = origin + stride*Bunch::y;
    ypa   = origin + stride*Bunch::yp;
    cdta  = origin + stride*Bunch::cdt;
    dpopa = origin + stride*Bunch::dpop;
}

#if 0
void Bunch::set_alt_arrays(double * RESTRICT &xa, double * RESTRICT &xpa,
                       double * RESTRICT &ya, double * RESTRICT &ypa,
                       double * RESTRICT &cdta, double * RESTRICT &dpopa)
{
    double *origin = alt_local_particles->origin();
    int stride = alt_local_particles->shape()[0];
    xa = origin + stride*Bunch::x;
    xpa = origin + stride*Bunch::xp;
    ya = origin + stride*Bunch::y;
    ypa = origin + stride*Bunch::yp;
    cdta = origin + stride*Bunch::cdt;
    dpopa = origin + stride*Bunch::dpop;
}
#endif

