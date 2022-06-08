/*
04/05/2016 Christoph Hasse
DISCLAIMER:

This code is not sufficiently tested yet and still under heavy development!
See *.cu file for more details
*/

#pragma once

#include <goofit/GlobalCudaDefines.h>
#include <goofit/Variable.h>
#include <goofit/detail/Complex.h>

#include <thrust/device_vector.h>
#include <thrust/iterator/constant_iterator.h>

namespace GooFit {

class ResonancePdf;
class Amplitude;
struct ParameterContainer;

template <typename E>
constexpr auto enum_to_underlying(E e) -> typename std::underlying_type<E>::type {
    return static_cast<typename std::underlying_type<E>::type>(e);
}

__host__ __device__ auto inDalitz(
    const fptype &m12, const fptype &m13, const fptype &bigM, const fptype &dm1, const fptype &dm2, const fptype &dm3)
    -> bool;

__device__ auto getResonanceAmplitude(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) -> fpcomplex;

__device__ void get4Vecs(fptype *Vecs,
                         const fptype &m12,
                         const fptype &m34,
                         const fptype &cos12,
                         const fptype &cos34,
                         const fptype &phi,
                         const fptype M,
                         const fptype m1,
                         const fptype m2,
                         const fptype m3,
                         const fptype m4);

__device__ auto getmass(const unsigned int &pair,
                        fptype &d1,
                        fptype &d2,
                        const fptype *vecs,
                        const fptype &m1,
                        const fptype &m2,
                        const fptype &m3,
                        const fptype &m4) -> fptype;

// in case of 3 particles the first two are the resonance.
enum DP4Pair {
    M_12 = 0,
    M_34,
    M_13,
    M_14,
    M_23,
    M_24,
    M_12_3,
    M_13_2,
    M_23_1,
    M_12_4,
    M_14_2,
    M_24_1,
    M_13_4,
    M_14_3,
    M_34_1,
    M_23_4,
    M_24_3,
    M_34_2
};

std::ostream &operator<<(std::ostream &out, const DP4Pair &obj);

enum DaughterPair { PAIR_12 = 0, PAIR_13, PAIR_23 };

const int resonanceSize = 4; // Number of parameters to describe one resonance.
// Why not make this a static const member of ResonancePdf? Because the 'static'
// keyword (and 'extern' as well) interacts badly with some nvcc versions when the
// variable is used in device code.

struct DecayInfo3 {
    fptype motherMass;
    fptype daug1Mass;
    fptype daug2Mass;
    fptype daug3Mass;
    fptype meson_radius;
    fptype mother_meson_radius;

    std::vector<ResonancePdf *> resonances;
};

struct DecayInfo3t : public DecayInfo3 {
    Variable _tau;
    Variable _xmixing;
    Variable _ymixing;
    Variable _deltax;
    Variable _deltay;

    DecayInfo3t(Variable _tau, Variable _xmixing, Variable _ymixing, Variable _deltax, Variable _deltay)
        : _tau(_tau)
        , _xmixing(_xmixing)
        , _ymixing(_ymixing)
        , _deltax(_deltax)
        , _deltay(_deltay) {}
};

struct DecayInfo4 {
    std::vector<fptype> particle_masses;
    fptype meson_radius;

    std::vector<Amplitude *> amplitudes;
    std::vector<Amplitude *> amplitudes_B;
};

struct DecayInfo4t : public DecayInfo4 {
    Variable _tau;
    Variable _xmixing;
    Variable _ymixing;
    Variable _SqWStoRSrate;

    DecayInfo4t(Variable _tau, Variable _xmixing, Variable _ymixing, Variable _SqWStoRSrate)
        : _tau(_tau)
        , _xmixing(_xmixing)
        , _ymixing(_ymixing)
        , _SqWStoRSrate(_SqWStoRSrate) {}
};

// Copied from strided_range thrust example by Nathan Bell.
// Iterator to move forward by a specified number of steps
// in each iteration.
template <typename Iterator>
class strided_range {
  public:
    typedef typename thrust::iterator_difference<Iterator>::type difference_type;

    struct stride_functor : public thrust::unary_function<difference_type, difference_type> {
        difference_type stride;

        stride_functor(difference_type stride)
            : stride(stride) {}

        __host__ __device__ auto operator()(const difference_type &i) const -> difference_type { return stride * i; }
    };
    typedef typename thrust::counting_iterator<difference_type> CountingIterator;
    typedef typename thrust::transform_iterator<stride_functor, CountingIterator> TransformIterator;
    typedef typename thrust::permutation_iterator<Iterator, TransformIterator> PermutationIterator;

    // type of the strided_range iterator
    typedef PermutationIterator iterator;

    // construct strided_range for the range [first,last)
    strided_range(Iterator first, Iterator last, difference_type stride)
        : first(first)
        , last(last)
        , stride(stride) {}

    auto begin() const -> iterator {
        return PermutationIterator(first, TransformIterator(CountingIterator(0), stride_functor(stride)));
    }

    auto end() const -> iterator { return begin() + ((last - first) + (stride - 1)) / stride; }

  protected:
    Iterator first;
    Iterator last;
    difference_type stride;
};

// From 3 body T

// thrust::tuple can't go down the read-only cache pipeline, so we are creating a structure for this.
typedef struct {
    // note: combining these into a single transaction (double2) should improve memory performance
    fptype ai_real;
    fptype ai_imag;
    fptype bi_real;
    fptype bi_imag;
} WaveHolder_s;

typedef thrust::tuple<fptype, fptype, fptype, fptype> WaveHolder;
typedef thrust::tuple<fptype, fptype, fptype, fptype, fptype, fptype> ThreeComplex;

} // namespace GooFit
