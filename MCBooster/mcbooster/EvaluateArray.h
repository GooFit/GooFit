/*
 * EvaluateArray.h
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
/*! \file EvaluateArray.h
 * Evaluate helper template function implementation.
 */

#ifndef EVALUATEARRAY_H_
#define EVALUATEARRAY_H_


#include <mcbooster/Config.h>
#include <mcbooster/Vector3R.h>
#include <mcbooster/Vector4R.h>
#include <mcbooster/GContainers.h>
#include <mcbooster/GTypes.h>
#include <mcbooster/functors/Calculate.h>

namespace MCBooster
{

/** Template functor for calculate an array of variables over a given set of particles.
 * Template functor for evaluate an arbitrary function object over the a set of particles stored
 * in the device. The function is supposed to evaluate at once many variables and the results are returned to the
 * __host__ via a given VariableSet_h. Datasets with up to nine particles can be handled.
 */

template<typename CUSTOMFUNC>
void EvaluateArray(const CUSTOMFUNC funcObj, ParticlesSet_d &pset,
		VariableSet_h &varset)
{

	if (pset.size() > 10 || pset.size() < 2)
	{
		cout
				<< "Can not Calculate(Eval) more than a nine-particle invariant mass."
				<< endl;
		return;
	}

	GInt_t arrayWidth = varset.size();
	GLong_t numberEvents = varset[0]->size();

	RealVector_d dev_array(arrayWidth * numberEvents);

	strided_range<RealVector_d::iterator> it_array(dev_array.begin(),
			dev_array.end(), arrayWidth);

	if (pset.size() == 2)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];

		thrust::for_each(
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.begin(), dev_v0->begin(),
								dev_v1->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.end(), dev_v0->end(),
								dev_v1->end())),
				Calculate3<CUSTOMFUNC>(funcObj));

	}

	if (pset.size() == 3)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];

		thrust::for_each(
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.begin(), dev_v0->begin(),
								dev_v1->begin(), dev_v2->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.end(), dev_v0->end(),
								dev_v1->end(), dev_v2->end())),
				Calculate3<CUSTOMFUNC>(funcObj));

	}

	if (pset.size() == 4)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];
		mc_device_vector<Vector4R> *dev_v3 = pset[3];

		thrust::for_each(
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.begin(), dev_v0->begin(),
								dev_v1->begin(), dev_v2->begin(),
								dev_v3->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.end(), dev_v0->end(),
								dev_v1->end(), dev_v2->end(), dev_v3->end())),
				Calculate3<CUSTOMFUNC>(funcObj));

	}

	if (pset.size() == 5)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];
		mc_device_vector<Vector4R> *dev_v3 = pset[3];
		mc_device_vector<Vector4R> *dev_v4 = pset[4];

		thrust::for_each(
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.begin(), dev_v0->begin(),
								dev_v1->begin(), dev_v2->begin(),
								dev_v3->begin(), dev_v4->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.end(), dev_v0->end(),
								dev_v1->end(), dev_v2->end(), dev_v3->end(),
								dev_v4->end())),
				Calculate3<CUSTOMFUNC>(funcObj));

	}

	if (pset.size() == 6)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];
		mc_device_vector<Vector4R> *dev_v3 = pset[3];
		mc_device_vector<Vector4R> *dev_v4 = pset[4];
		mc_device_vector<Vector4R> *dev_v5 = pset[5];

		thrust::for_each(
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.begin(), dev_v0->begin(),
								dev_v1->begin(), dev_v2->begin(),
								dev_v3->begin(), dev_v4->begin(),
								dev_v5->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.end(), dev_v0->end(),
								dev_v1->end(), dev_v2->end(), dev_v3->end(),
								dev_v4->end(), dev_v5->end())),
				Calculate3<CUSTOMFUNC>(funcObj));

	}

	if (pset.size() == 7)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];
		mc_device_vector<Vector4R> *dev_v3 = pset[3];
		mc_device_vector<Vector4R> *dev_v4 = pset[4];
		mc_device_vector<Vector4R> *dev_v5 = pset[5];
		mc_device_vector<Vector4R> *dev_v6 = pset[6];

		thrust::for_each(
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.begin(), dev_v0->begin(),
								dev_v1->begin(), dev_v2->begin(),
								dev_v3->begin(), dev_v4->begin(),
								dev_v5->begin(), dev_v6->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.end(), dev_v0->end(),
								dev_v1->end(), dev_v2->end(), dev_v3->end(),
								dev_v4->end(), dev_v5->end(), dev_v6->end())),
				Calculate3<CUSTOMFUNC>(funcObj));

	}

	if (pset.size() == 8)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];
		mc_device_vector<Vector4R> *dev_v3 = pset[3];
		mc_device_vector<Vector4R> *dev_v4 = pset[4];
		mc_device_vector<Vector4R> *dev_v5 = pset[5];
		mc_device_vector<Vector4R> *dev_v6 = pset[6];
		mc_device_vector<Vector4R> *dev_v7 = pset[7];

		thrust::for_each(
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.begin(), dev_v0->begin(),
								dev_v1->begin(), dev_v2->begin(),
								dev_v3->begin(), dev_v4->begin(),
								dev_v5->begin(), dev_v6->begin(),
								dev_v7->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.end(), dev_v0->end(),
								dev_v1->end(), dev_v2->end(), dev_v3->end(),
								dev_v4->end(), dev_v5->end(), dev_v6->end(),
								dev_v7->end())),
				Calculate3<CUSTOMFUNC>(funcObj));

	}

	if (pset.size() == 9)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];
		mc_device_vector<Vector4R> *dev_v3 = pset[3];
		mc_device_vector<Vector4R> *dev_v4 = pset[4];
		mc_device_vector<Vector4R> *dev_v5 = pset[5];
		mc_device_vector<Vector4R> *dev_v6 = pset[6];
		mc_device_vector<Vector4R> *dev_v7 = pset[7];
		mc_device_vector<Vector4R> *dev_v8 = pset[8];

		thrust::for_each(
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.begin(), dev_v0->begin(),
								dev_v1->begin(), dev_v2->begin(),
								dev_v3->begin(), dev_v4->begin(),
								dev_v5->begin(), dev_v6->begin(),
								dev_v7->begin(), dev_v8->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.end(), dev_v0->end(),
								dev_v1->end(), dev_v2->end(), dev_v3->end(),
								dev_v4->end(), dev_v5->end(), dev_v6->end(),
								dev_v7->end(), dev_v8->end())),
				Calculate3<CUSTOMFUNC>(funcObj));

	}

	if (pset.size() == 10)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];
		mc_device_vector<Vector4R> *dev_v3 = pset[3];
		mc_device_vector<Vector4R> *dev_v4 = pset[4];
		mc_device_vector<Vector4R> *dev_v5 = pset[5];
		mc_device_vector<Vector4R> *dev_v6 = pset[6];
		mc_device_vector<Vector4R> *dev_v7 = pset[7];
		mc_device_vector<Vector4R> *dev_v8 = pset[8];

		thrust::for_each(
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.begin(), dev_v0->begin(),
								dev_v1->begin(), dev_v2->begin(),
								dev_v3->begin(), dev_v4->begin(),
								dev_v5->begin(), dev_v6->begin(),
								dev_v7->begin(), dev_v8->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.end(), dev_v0->end(),
								dev_v1->end(), dev_v2->end(), dev_v3->end(),
								dev_v4->end(), dev_v5->end(), dev_v6->end(),
								dev_v7->end(), dev_v8->end())),
				Calculate3<CUSTOMFUNC>(funcObj));

	}

#if THRUST_DEVICE_SYSTEM==THRUST_DEVICE_BACKEND_OMP || THRUST_DEVICE_SYSTEM==THRUST_DEVICE_BACKEND_TBB

#pragma omp parallel num_threads(  arrayWidth )
	{
		strided_range<RealVector_d::iterator> it_array(dev_array.begin() + omp_get_thread_num()
				, dev_array.end(), arrayWidth);

		thrust::copy(it_array.begin(),it_array.end(),
				varset[omp_get_thread_num()]->begin());
	}
#else
	cudaStream_t s[arrayWidth];

	for (GInt_t d = 0; d < arrayWidth; d++)
	{
		cudaStreamCreate(&s[d]);

	}
	strided_range<RealVector_d::iterator> *it[arrayWidth];
	for (GInt_t d = 0; d < arrayWidth; d++)
		it[d] = new strided_range<RealVector_d::iterator>(dev_array.begin() + d,
				dev_array.end(), arrayWidth);
	for (GInt_t d = 0; d < arrayWidth; d++)
	{

		thrust::copy(thrust::cuda::par.on(s[d]), it[d]->begin(), it[d]->end(),
				varset[d]->begin());

	}
	cudaDeviceSynchronize();
	for (GInt_t d = 0; d < arrayWidth; d++)
		cudaStreamDestroy(s[d]);
	for (GInt_t d = 0; d < arrayWidth; d++)
		delete it[d];

#endif
	return;
}

#if !(THRUST_DEVICE_SYSTEM==THRUST_DEVICE_BACKEND_OMP || THRUST_DEVICE_SYSTEM==THRUST_DEVICE_BACKEND_TBB)

/** Template functor for calculate an array of variables over a given set of particles.
 * Template functor for evaluate an arbitrary function object over the a set of particles stored
 * in the device. The function is supposed to evaluate at once many variables and the results are returned to the
 * __device__ via a given VariableSet_h. Datasets with up to nine particles  can be handled. __This function is available
 * only for CUDA backends.__
 */

template<typename CUSTOMFUNC>
void EvaluateArray(const CUSTOMFUNC funcObj, ParticlesSet_d &pset,
		VariableSet_d &varset)
{

	if (pset.size() > 10 || pset.size() < 2)
	{
		cout
				<< "Can not Calculate(Eval) more than a nine-particle invariant mass."
				<< endl;
		return;
	}

	GInt_t arrayWidth = varset.size();
	GLong_t numberEvents = varset[0]->size();

	RealVector_d dev_array(arrayWidth * numberEvents);

	strided_range<RealVector_d::iterator> it_array(dev_array.begin(),
			dev_array.end(), arrayWidth);

	if (pset.size() == 2)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];

		thrust::for_each(
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.begin(), dev_v0->begin(),
								dev_v1->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.end(), dev_v0->end(),
								dev_v1->end())),
				Calculate3<CUSTOMFUNC>(funcObj));

	}

	if (pset.size() == 3)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];

		thrust::for_each(
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.begin(), dev_v0->begin(),
								dev_v1->begin(), dev_v2->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.end(), dev_v0->end(),
								dev_v1->end(), dev_v2->end())),
				Calculate3<CUSTOMFUNC>(funcObj));

	}

	if (pset.size() == 4)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];
		mc_device_vector<Vector4R> *dev_v3 = pset[3];

		thrust::for_each(
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.begin(), dev_v0->begin(),
								dev_v1->begin(), dev_v2->begin(),
								dev_v3->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.end(), dev_v0->end(),
								dev_v1->end(), dev_v2->end(), dev_v3->end())),
				Calculate3<CUSTOMFUNC>(funcObj));

	}

	if (pset.size() == 5)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];
		mc_device_vector<Vector4R> *dev_v3 = pset[3];
		mc_device_vector<Vector4R> *dev_v4 = pset[4];

		thrust::for_each(
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.begin(), dev_v0->begin(),
								dev_v1->begin(), dev_v2->begin(),
								dev_v3->begin(), dev_v4->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.end(), dev_v0->end(),
								dev_v1->end(), dev_v2->end(), dev_v3->end(),
								dev_v4->end())),
				Calculate3<CUSTOMFUNC>(funcObj));

	}

	if (pset.size() == 6)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];
		mc_device_vector<Vector4R> *dev_v3 = pset[3];
		mc_device_vector<Vector4R> *dev_v4 = pset[4];
		mc_device_vector<Vector4R> *dev_v5 = pset[5];

		thrust::for_each(
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.begin(), dev_v0->begin(),
								dev_v1->begin(), dev_v2->begin(),
								dev_v3->begin(), dev_v4->begin(),
								dev_v5->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.end(), dev_v0->end(),
								dev_v1->end(), dev_v2->end(), dev_v3->end(),
								dev_v4->end(), dev_v5->end())),
				Calculate3<CUSTOMFUNC>(funcObj));

	}

	if (pset.size() == 7)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];
		mc_device_vector<Vector4R> *dev_v3 = pset[3];
		mc_device_vector<Vector4R> *dev_v4 = pset[4];
		mc_device_vector<Vector4R> *dev_v5 = pset[5];
		mc_device_vector<Vector4R> *dev_v6 = pset[6];

		thrust::for_each(
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.begin(), dev_v0->begin(),
								dev_v1->begin(), dev_v2->begin(),
								dev_v3->begin(), dev_v4->begin(),
								dev_v5->begin(), dev_v6->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.end(), dev_v0->end(),
								dev_v1->end(), dev_v2->end(), dev_v3->end(),
								dev_v4->end(), dev_v5->end(), dev_v6->end())),
				Calculate3<CUSTOMFUNC>(funcObj));

	}

	if (pset.size() == 8)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];
		mc_device_vector<Vector4R> *dev_v3 = pset[3];
		mc_device_vector<Vector4R> *dev_v4 = pset[4];
		mc_device_vector<Vector4R> *dev_v5 = pset[5];
		mc_device_vector<Vector4R> *dev_v6 = pset[6];
		mc_device_vector<Vector4R> *dev_v7 = pset[7];

		thrust::for_each(
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.begin(), dev_v0->begin(),
								dev_v1->begin(), dev_v2->begin(),
								dev_v3->begin(), dev_v4->begin(),
								dev_v5->begin(), dev_v6->begin(),
								dev_v7->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.end(), dev_v0->end(),
								dev_v1->end(), dev_v2->end(), dev_v3->end(),
								dev_v4->end(), dev_v5->end(), dev_v6->end(),
								dev_v7->end())),
				Calculate3<CUSTOMFUNC>(funcObj));

	}

	if (pset.size() == 9)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];
		mc_device_vector<Vector4R> *dev_v3 = pset[3];
		mc_device_vector<Vector4R> *dev_v4 = pset[4];
		mc_device_vector<Vector4R> *dev_v5 = pset[5];
		mc_device_vector<Vector4R> *dev_v6 = pset[6];
		mc_device_vector<Vector4R> *dev_v7 = pset[7];
		mc_device_vector<Vector4R> *dev_v8 = pset[8];

		thrust::for_each(
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.begin(), dev_v0->begin(),
								dev_v1->begin(), dev_v2->begin(),
								dev_v3->begin(), dev_v4->begin(),
								dev_v5->begin(), dev_v6->begin(),
								dev_v7->begin(), dev_v8->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.end(), dev_v0->end(),
								dev_v1->end(), dev_v2->end(), dev_v3->end(),
								dev_v4->end(), dev_v5->end(), dev_v6->end(),
								dev_v7->end(), dev_v8->end())),
				Calculate3<CUSTOMFUNC>(funcObj));

	}

	if (pset.size() == 10)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];
		mc_device_vector<Vector4R> *dev_v3 = pset[3];
		mc_device_vector<Vector4R> *dev_v4 = pset[4];
		mc_device_vector<Vector4R> *dev_v5 = pset[5];
		mc_device_vector<Vector4R> *dev_v6 = pset[6];
		mc_device_vector<Vector4R> *dev_v7 = pset[7];
		mc_device_vector<Vector4R> *dev_v8 = pset[8];

		thrust::for_each(
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.begin(), dev_v0->begin(),
								dev_v1->begin(), dev_v2->begin(),
								dev_v3->begin(), dev_v4->begin(),
								dev_v5->begin(), dev_v6->begin(),
								dev_v7->begin(), dev_v8->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(it_array.end(), dev_v0->end(),
								dev_v1->end(), dev_v2->end(), dev_v3->end(),
								dev_v4->end(), dev_v5->end(), dev_v6->end(),
								dev_v7->end(), dev_v8->end())),
				Calculate3<CUSTOMFUNC>(funcObj));

	}

#if THRUST_DEVICE_SYSTEM==THRUST_DEVICE_BACKEND_OMP || THRUST_DEVICE_SYSTEM==THRUST_DEVICE_BACKEND_TBB

#pragma omp parallel num_threads(  arrayWidth )
	{
		strided_range<RealVector_d::iterator> it_array(dev_array.begin() + omp_get_thread_num()
				, dev_array.end(), arrayWidth);

		thrust::copy(it_array.begin(),it_array.end(),
				varset[omp_get_thread_num()]->begin());
	}
#else
	cudaStream_t s[arrayWidth];

	for (GInt_t d = 0; d < arrayWidth; d++)
	{
		cudaStreamCreate(&s[d]);

	}
	strided_range<RealVector_d::iterator> *it[arrayWidth];
	for (GInt_t d = 0; d < arrayWidth; d++)
		it[d] = new strided_range<RealVector_d::iterator>(dev_array.begin() + d,
				dev_array.end(), arrayWidth);
	for (GInt_t d = 0; d < arrayWidth; d++)
	{

		thrust::copy(thrust::cuda::par.on(s[d]), it[d]->begin(), it[d]->end(),
				varset[d]->begin());

	}
	cudaDeviceSynchronize();
	for (GInt_t d = 0; d < arrayWidth; d++)
		cudaStreamDestroy(s[d]);
	for (GInt_t d = 0; d < arrayWidth; d++)
		delete it[d];

#endif
	return;
}
#endif

}

#endif /* EVALUATEARRAY_H_ */
