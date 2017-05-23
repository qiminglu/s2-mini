#ifndef MULTI_ARRAY_TYPEDEFS_H_
#define MULTI_ARRAY_TYPEDEFS_H_

#include <complex>
#include <vector>


// Synergia Vector (svec)
template<typename T>
class svec
{
public:

    svec()
      : data(new T[0]), size(0)
    {
        #pragma acc enter data copyin(this[0:1]) create(data[0:size])
    }

    svec(size_t size_)
      : data(new T[size_]), size(size_)
    {
        #pragma acc enter data copyin(this[0:1]) create(data[0:size])
    }

    ~svec()
    {
        #pragma acc exit data delete(data[0:size], this[0:1])
        delete [] data;
    }

    void resize(size_t size_)
    {
        delete [] data;
        data = new T[size_];
        size = size_;

        #pragma acc enter data copyin(this[0:1]) create(data[0:size])
    }

    void copy_host_to_dev()
    {
        #pragma acc update device(data[0:size])
    }

    void copy_dev_to_host()
    {
        #pragma acc update self(data[0:size])
    }

    inline T & operator[] (size_t i) const
    { return data[i]; }

    inline T * origin() const
    { return data; }

private:

    T    * data;
    size_t size;
};


typedef std::vector<double>         MArray1d;
typedef std::vector<double> &       MArray1d_ref;
typedef std::vector<double> const & Const_MArray1d_ref;

#if 0
typedef std::vector<double>         MArray2d;
typedef std::vector<double> &       MArray2d_ref;
typedef std::vector<double> const & Const_MArray2d_ref;
#endif

typedef svec<double>                MArray2d;
typedef svec<double> &              MArray2d_ref;
typedef svec<double> const &        Const_MArray2d_ref;


#if 0
#include "boost/multi_array.hpp"
#include "boost/shared_array.hpp"

#include "boost/align/aligned_allocator.hpp"
#endif

#if 0
typedef boost::multi_array_types::index_range range;
typedef boost::multi_array_types::extent_range extent_range;

template<typename T, size_t N_dims>
    struct Raw_multi_array
    {
        boost::const_multi_array_ref<T, N_dims > * dummy;
        typedef typename boost::const_multi_array_ref<T, N_dims >::storage_order_type storage_order_type;
        boost::shared_array<T > p;
        boost::multi_array_ref<T, N_dims > m;
        template<class ExtentList>
            explicit
            Raw_multi_array(const ExtentList& extents,
                    const storage_order_type& store = boost::c_storage_order()) :
                    dummy(
                            new boost::const_multi_array_ref<T, N_dims >(
                                    static_cast<T* >(0L), extents, store)), p(
                            new T[dummy->num_elements()]), m(p.get(), extents, store)
            {
                delete dummy;
            }
        ~Raw_multi_array()
        {
        }
    };

typedef boost::multi_array<double, 1 > MArray1d; // syndoc:include
typedef boost::multi_array_ref<double, 1 > MArray1d_ref; // syndoc:include
typedef boost::const_multi_array_ref<double, 1 > Const_MArray1d_ref; // syndoc:include
typedef boost::detail::multi_array::multi_array_view<double, 1 > MArray1d_view;  // syndoc:include
typedef Raw_multi_array<double, 1> Raw_MArray1d; // syndoc:include


typedef boost::multi_array<double, 2 > MArray2d; // syndoc:include
typedef boost::multi_array<double, 2, boost::alignment::aligned_allocator<double, 64> > MArray2da; // syndoc:include
typedef boost::multi_array_ref<double, 2 > MArray2d_ref; // syndoc:include
typedef boost::const_multi_array_ref<double, 2 > Const_MArray2d_ref; // syndoc:include
typedef boost::detail::multi_array::multi_array_view<double, 2 > MArray2d_view; // syndoc:include
typedef boost::detail::multi_array::const_multi_array_view<double, 2 >
        Const_MArray2d_view; // syndoc:include
typedef boost::general_storage_order<2> storage2d;
typedef Raw_multi_array<double, 2> Raw_MArray2d; // syndoc:include

typedef boost::multi_array<double, 3 > MArray3d; // syndoc:include
typedef boost::multi_array_ref<double, 3 > MArray3d_ref; // syndoc:include
typedef boost::const_multi_array_ref<double, 3 > Const_MArray3d_ref; // syndoc:include
typedef boost::detail::multi_array::multi_array_view<double, 3 > MArray3d_view; // syndoc:include
typedef boost::general_storage_order<3> storage3d;
typedef Raw_multi_array<double, 3> Raw_MArray3d; // syndoc:include

typedef boost::multi_array<std::complex<double >, 1 > MArray1dc; // syndoc:include
typedef boost::multi_array_ref<std::complex<double >, 1 > MArray1dc_ref; // syndoc:include
typedef boost::const_multi_array_ref<std::complex<double >, 1 >
        Const_MArray1dc_ref; // syndoc:include
typedef boost::detail::multi_array::multi_array_view<std::complex<double >, 1 >
        MArray1dc_view; // syndoc:include

typedef boost::multi_array<std::complex<double >, 2 > MArray2dc; // syndoc:include
typedef boost::multi_array_ref<std::complex<double >, 2 > MArray2dc_ref; // syndoc:include
typedef boost::const_multi_array_ref<std::complex<double >, 2 >
        Const_MArray2dc_ref; // syndoc:include
typedef boost::detail::multi_array::multi_array_view<std::complex<double >, 2 >
        MArray2dc_view; // syndoc:include
typedef Raw_multi_array<std::complex<double >, 2> Raw_MArray2dc; // syndoc:include

typedef boost::multi_array<std::complex<double >, 3 > MArray3dc; // syndoc:include
typedef boost::multi_array_ref<std::complex<double >, 3 > MArray3dc_ref; // syndoc:include
typedef boost::const_multi_array_ref<std::complex<double >, 3 >
        Const_MArray3dc_ref; // syndoc:include
typedef boost::detail::multi_array::multi_array_view<std::complex<double >, 3 >
        MArray3dc_view; // syndoc:include

typedef boost::multi_array<int, 1 > MArray1i; // syndoc:include
typedef boost::multi_array_ref<int, 1 > MArray1i_ref; // syndoc:include
typedef boost::const_multi_array_ref<int, 1 > Const_MArray1i_ref; // syndoc:include
typedef boost::detail::multi_array::multi_array_view<int, 1 > MArray1i_view; // syndoc:include
#endif

#endif /* MULTI_ARRAY_TYPEDEFS_H_ */
