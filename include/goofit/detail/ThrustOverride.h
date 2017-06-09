#pragma once

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA

#include <thrust/execution_policy.h>
#include <thrust/detail/config.h>
#include <thrust/reduce.h>
#include <thrust/detail/seq.h>
#include <thrust/detail/temporary_array.h>
#include <thrust/system/cuda/detail/bulk.h>
#include <thrust/system/cuda/detail/decomposition.h>
#include <thrust/system/cuda/detail/execution_policy.h>
#include <thrust/system/cuda/detail/execute_on_stream.h>
#include <thrust/detail/type_traits.h>

#include <thrust/system/cuda/detail/bulk.h>

#include "goofit/detail/ThrustOverrideConfig.h"

#ifndef OR_GROUPSIZE
#define OR_GROUPSIZE 128
#endif

#ifndef OR_GRAINSIZE
#define OR_GRAINSIZE 7
#endif

// We are defining a new policy.  This will allow us to override the internal thread/block distribution
struct goofit_policy : thrust::device_execution_policy<goofit_policy> {};

template<typename InputIterator, typename UnaryFunction, typename OutputType, typename BinaryFunction>
__host__ __device__ OutputType transform_reduce(goofit_policy &exec,
                                                InputIterator first,
                                                InputIterator last,
                                                UnaryFunction unary_op,
                                                OutputType init,
                                                BinaryFunction binary_op) {
    // printf ("goofit transform_reduce\n");
    thrust::transform_iterator<UnaryFunction, InputIterator, OutputType> xfrm_first(first, unary_op);
    thrust::transform_iterator<UnaryFunction, InputIterator, OutputType> xfrm_last(last, unary_op);

    return thrust::reduce(exec, xfrm_first, xfrm_last, init, binary_op);
} // end transform_reduce()

template<typename InputIterator, typename OutputType, typename BinaryFunction>
__host__ __device__ OutputType
reduce(goofit_policy &exec, InputIterator first, InputIterator last, OutputType init, BinaryFunction binary_op) {
    typedef typename thrust::iterator_difference<InputIterator>::type size_type;

    const size_type n = last - first;

    if(n <= 0)
        return init;

    // We are no longer expecting a default stream to be passed
    cudaStream_t s = stream(thrust::detail::derived_cast(exec));

    // OR_GROUPSIZE=128
    // OR_GRAINSIZE=7
    // printf ("groupsize:%i grainsize:%i\n", OR_GROUPSIZE, OR_GRAINSIZE);
    const size_type groupsize = OR_GROUPSIZE;
    // const size_type groupsize = 128;
    // const size_type grainsize = 7;
    const size_type grainsize    = OR_GRAINSIZE;
    const size_type tile_size    = groupsize * grainsize;
    const size_type num_tiles    = (n + tile_size - 1) / tile_size;
    const size_type subscription = 10;

    thrust::system::cuda::detail::bulk_::concurrent_group<thrust::system::cuda::detail::bulk_::agent<grainsize>,
                                                          groupsize>
        g;

    const size_type num_groups = thrust::min<size_type>(subscription * g.hardware_concurrency(), num_tiles);

    thrust::system::cuda::detail::aligned_decomposition<size_type> decomp(n, num_groups, tile_size);

    thrust::detail::temporary_array<OutputType, goofit_policy> partial_sums(exec, decomp.size());

    // reduce into partial sums
    thrust::system::cuda::detail::bulk_::async(thrust::system::cuda::detail::bulk_::par(s, g, decomp.size()),
                                               thrust::system::cuda::detail::reduce_detail::reduce_partitions(),
                                               thrust::system::cuda::detail::bulk_::root.this_exec,
                                               first,
                                               decomp,
                                               partial_sums.begin(),
                                               init,
                                               binary_op)
        .wait();

    if(partial_sums.size() > 1) {
        // reduce the partial sums
        thrust::system::cuda::detail::bulk_::async(thrust::system::cuda::detail::bulk_::par(s, g, 1),
                                                   thrust::system::cuda::detail::reduce_detail::reduce_partitions(),
                                                   thrust::system::cuda::detail::bulk_::root.this_exec,
                                                   partial_sums.begin(),
                                                   partial_sums.end(),
                                                   partial_sums.begin(),
                                                   binary_op);
    } // end if

    return get_value(exec, &partial_sums[0]);
} // end goofit_reduce()

#endif
