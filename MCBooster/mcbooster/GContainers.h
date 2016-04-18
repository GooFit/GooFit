/*
 * GContainers.h
 *
 * Copyright 2016 Antonio Augusto Alves Junior
 *
 * Created on : Feb 25, 2016
 *      Author: Antonio Augusto Alves Junior
 */

/*
 *   This file is part of MCBooster.
 *
 *   MCBooster is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   MCBooster is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with MCBooster.  If not, see <http://www.gnu.org/licenses/>.
 */

/*! \file GContainers.h
 *  Typedef for useful container classes used in MCBooster.
 *  Containers defined here should be used in users application also.
 */

#ifndef GCONTAINERS_H_
#define GCONTAINERS_H_


#include <mcbooster/Config.h>
#include <mcbooster/Vector3R.h>
#include <mcbooster/Vector4R.h>
#include <vector>
#include <mcbooster/GTypes.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/complex.h>

#if !(THRUST_DEVICE_SYSTEM==THRUST_DEVICE_BACKEND_OMP || THRUST_DEVICE_SYSTEM==THRUST_DEVICE_BACKEND_TBB)
#include <thrust/system/cuda/experimental/pinned_allocator.h>
#endif

using namespace std;

namespace MCBooster
{

#if (THRUST_DEVICE_SYSTEM==THRUST_DEVICE_BACKEND_OMP || THRUST_DEVICE_SYSTEM==THRUST_DEVICE_BACKEND_TBB)
/*!
 * Generic template typedef for thrust::host_vector. Use it instead of Thrust implementation
 * in order to avoid problems to compile OpenMP based applications using gcc and without a cuda runtime installation.
 */
	template <typename T>
		using  mc_device_vector = thrust::host_vector<T>;
/*!
 * Generic template typedef for thrust::host_vector. Use it instead of Thrust implementation
 * in order to avoid problems to compile OpenMP based applications using gcc and without a cuda runtime installation.
 * mc_host_vectot will always allocate page locked memory on CUDA backends in order to maximize speed in memory transfers
 * to the device.
 */
	template <typename T>
		using  mc_host_vector   = thrust::host_vector<T>;

#else
	/*!
	 * Generic template typedef for thrust::host_vector. Use it instead of Thrust implementation
	 * in order to avoid problems to compile OpenMP based applications using gcc and without a cuda runtime installation.
	 */
	template <typename T>
		using  mc_device_vector = thrust::device_vector<T>;
	/*!
	 * Generic template typedef for thrust::host_vector. Use it instead of Thrust implementation
	 * in order to avoid problems to compile OpenMP based applications using gcc and without a cuda runtime installation.
	 * mc_host_vectot will always allocate page locked memory on CUDA backends in order to maximize speed in memory transfers
	 * to the device.
	 */
	template <typename T>
		using  mc_host_vector = thrust::host_vector<T,
				thrust::cuda::experimental::pinned_allocator<T>>;

#endif


//-----------------------------------------------------------------------
//complex number container
typedef thrust::complex<GReal_t> GComplex_t; /*! Typedef for complex number.*/

//-----------------------------------------------------------------------

typedef mc_host_vector<Vector4R> FourVectors_h; /*! Vector4R host vector. Use it to store four-vectors at __host__.*/
typedef mc_host_vector<Vector3R> ThreeVectors_h; /*! Vector3R host vector. Use it to store four-vectors at __host__.*/

//-----------------------------------------------------------------------
//basic containers on host

typedef mc_host_vector<GBool_t>    BoolVector_h; /*! Typedef for a GBool_t host vector.*/
typedef mc_host_vector<GReal_t>    RealVector_h; /*! Typedef for a GReal_t host vector.*/
typedef mc_host_vector<GComplex_t> ComplexVector_h; /*! Typedef for a GComplex_t host vector.*/
typedef mc_host_vector<Vector4R>   Particles_h;/*! Typedef for a  Vector4R host vector.*/
typedef vector<Particles_h*>       ParticlesSet_h;/*! Typedef for a  STL vector of pointers to host Particles_h vectors .*/
typedef vector<RealVector_h*>      VariableSet_h;/*! Typedef for a STL vector of pointers to host RealVector_h vectors.*/


//-----------------------------------------------------------------------
//basic containers on device
typedef mc_device_vector<GBool_t>    BoolVector_d; /*! Typedef for a GBool_t device vector.*/
typedef mc_device_vector<GReal_t>    RealVector_d; /*! Typedef for a GReal_t device vector.*/
typedef mc_device_vector<GComplex_t> ComplexVector_d; /*! Typedef for a GComplex_t device vector.*/
typedef mc_device_vector<Vector4R>   Particles_d; /*! Typedef for a  Vector4R device vector.*/
typedef vector<Particles_d*>  ParticlesSet_d; /*! Typedef for a  STL vector of pointers to device Particles_d vectors.*/
typedef vector<RealVector_d*> VariableSet_d; /*! Typedef for a STL vector of pointers to device RealVector_d vectors.*/


/*! GT1 iterator is a typedef for thrust::detail::tuple_of_iterator_references<Vector4R &,...> */
typedef thrust::detail::tuple_of_iterator_references<Vector4R &,
		thrust::null_type, thrust::null_type, thrust::null_type,
		thrust::null_type, thrust::null_type, thrust::null_type,
		thrust::null_type, thrust::null_type, thrust::null_type> GT1;

/*! GT2 iterator is a typedef for thrust::detail::tuple_of_iterator_references<Vector4R &, Vector4R &, ...> */
typedef thrust::detail::tuple_of_iterator_references<Vector4R &, Vector4R &,
		thrust::null_type, thrust::null_type, thrust::null_type,
		thrust::null_type, thrust::null_type, thrust::null_type,
		thrust::null_type, thrust::null_type> GT2;

/*! GT3 iterator is a typedef for thrust::detail::tuple_of_iterator_references<Vector4R &, Vector4R &, Vector4R &, ...> */
typedef thrust::detail::tuple_of_iterator_references<Vector4R &, Vector4R &,
		Vector4R &, thrust::null_type, thrust::null_type, thrust::null_type,
		thrust::null_type, thrust::null_type, thrust::null_type,
		thrust::null_type> GT3;

/*! GT4  iterator is a typedef for thrust::detail::tuple_of_iterator_references<Vector4R &, Vector4R &, Vector4R &, Vector4R &,...> */
typedef thrust::detail::tuple_of_iterator_references<Vector4R &, Vector4R &,
		Vector4R &, Vector4R &, thrust::null_type, thrust::null_type,
		thrust::null_type, thrust::null_type, thrust::null_type,
		thrust::null_type> GT4;

/*! GT5 iterator is a typedef for thrust::detail::tuple_of_iterator_references<Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &,...>*/
typedef thrust::detail::tuple_of_iterator_references<Vector4R &, Vector4R &,
		Vector4R &, Vector4R &, Vector4R &, thrust::null_type,
		thrust::null_type, thrust::null_type, thrust::null_type,
		thrust::null_type> GT5;

/*! GT6 iterator is a typedef for thrust::detail::tuple_of_iterator_references<Vector4R &, Vector4R &, Vector4R &, Vector4R &,
  Vector4R &,Vector4R &,...> */
typedef thrust::detail::tuple_of_iterator_references<Vector4R &, Vector4R &,
		Vector4R &, Vector4R &, Vector4R &, Vector4R &, thrust::null_type,
		thrust::null_type, thrust::null_type, thrust::null_type> GT6;

/*! GT7  iterator is a typedef for thrust::detail::tuple_of_iterator_references<Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &,Vector4R &,Vector4R &,...> */
typedef thrust::detail::tuple_of_iterator_references<Vector4R &, Vector4R &,
		Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &,
		thrust::null_type, thrust::null_type, thrust::null_type> GT7;

/*!GT8 iterator is a typedef for thrust::detail::tuple_of_iterator_references<Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &,Vector4R &,Vector4R &,Vector4R &,...> */
typedef thrust::detail::tuple_of_iterator_references<Vector4R &, Vector4R &,
		Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &,
		thrust::null_type, thrust::null_type> GT8;

/*! GT9  iterator is a typedef for thrust::detail::tuple_of_iterator_references<Vector4R &, Vector4R &, Vector4R &, Vector4R &,
   Vector4R &,Vector4R &,Vector4R &,Vector4R &,Vector4R &...> */
typedef thrust::detail::tuple_of_iterator_references<Vector4R &, Vector4R &,
		Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &,
		Vector4R &, thrust::null_type> GT9;

/*! GT10 iterator is a typedef for thrust::detail::tuple_of_iterator_references<Vector4R &, Vector4R &, Vector4R &, Vector4R &,  Vector4R &,Vector4R &,Vector4R &,Vector4R &,Vector4R &,Vector4R &> */
typedef thrust::detail::tuple_of_iterator_references<Vector4R &, Vector4R &,
		Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &,
		Vector4R &, Vector4R &> GT10;

/*!
 *
 * GTR2  iterator is a typedef for thrust::detail::tuple_of_iterator_references<GReal_t &, Vector4R &, ...>
 */
typedef thrust::detail::tuple_of_iterator_references<GReal_t &, Vector4R &,
		thrust::null_type, thrust::null_type, thrust::null_type,
		thrust::null_type, thrust::null_type, thrust::null_type,
		thrust::null_type, thrust::null_type> GTR2;

/*!
 *
 * GTR3 iterator is a typedef for thrust::detail::tuple_of_iterator_references<GReal_t &, Vector4R &, Vector4R &, ...>
 */
typedef thrust::detail::tuple_of_iterator_references<GReal_t &, Vector4R &,
		Vector4R &, thrust::null_type, thrust::null_type, thrust::null_type,
		thrust::null_type, thrust::null_type, thrust::null_type,
		thrust::null_type> GTR3;

/*!
 *
 * GTR4 iterator is a typedef for thrust::detail::tuple_of_iterator_references<GReal_t &, Vector4R &, Vector4R &, Vector4R &, ...>
 */
typedef thrust::detail::tuple_of_iterator_references<GReal_t &, Vector4R &,
		Vector4R &, Vector4R &, thrust::null_type, thrust::null_type,
		thrust::null_type, thrust::null_type, thrust::null_type,
		thrust::null_type> GTR4;
/*!
 *
 * GTR5 iterator is a typedef for thrust::detail::tuple_of_iterator_references<GReal_t &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, ...>
 */
typedef thrust::detail::tuple_of_iterator_references<GReal_t &, Vector4R &,
		Vector4R &, Vector4R &, Vector4R &, thrust::null_type,
		thrust::null_type, thrust::null_type, thrust::null_type,
		thrust::null_type> GTR5;

/*! GTR6 iterator is a typedef for thrust::detail::tuple_of_iterator_references<GReal_t &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, ...>*/
typedef thrust::detail::tuple_of_iterator_references<GReal_t &, Vector4R &,
		Vector4R &, Vector4R &, Vector4R &, Vector4R &, thrust::null_type,
		thrust::null_type, thrust::null_type, thrust::null_type> GTR6;

/*! GTR7 iterator is a typedef for thrust::detail::tuple_of_iterator_references<GReal_t &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, ...>*/
typedef thrust::detail::tuple_of_iterator_references<GReal_t &, Vector4R &,
		Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &,
		thrust::null_type, thrust::null_type, thrust::null_type> GTR7;

/*! GTR8 iterator is a typedef for thrust::detail::tuple_of_iterator_references<GReal_t &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, ...>*/
typedef thrust::detail::tuple_of_iterator_references<GReal_t &, Vector4R &,
		Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &,
		thrust::null_type, thrust::null_type> GTR8;

/*!GTR9  iterator is a typedef for thrust::detail::tuple_of_iterator_references<GReal_t &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, ...>*/
typedef thrust::detail::tuple_of_iterator_references<GReal_t &, Vector4R &,
		Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &,
		Vector4R &, thrust::null_type> GTR9;

/*! GTR10 iterator is a typedef for thrust::detail::tuple_of_iterator_references<GReal_t &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, ...>*/
typedef thrust::detail::tuple_of_iterator_references<GReal_t &, Vector4R &,
		Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &, Vector4R &,
		Vector4R &, Vector4R &> GTR10;

}
#endif /* GCONTAINERS_H_ */
