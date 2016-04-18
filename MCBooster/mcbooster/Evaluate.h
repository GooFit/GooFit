/*
 * Evaluate.h
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

/*! \file Evaluate.h
 * Evaluate helper template function implementation.
 */

#ifndef EVALUATE_H_
#define EVALUATE_H_


#include <mcbooster/Config.h>
#include <mcbooster/Vector3R.h>
#include <mcbooster/Vector4R.h>
#include <mcbooster/GContainers.h>
#include <mcbooster/GTypes.h>
#include <mcbooster/functors/Calculate.h>

namespace MCBooster
{
/** Template functor for evaluate an arbitrary function object.
 * Template functor for evaluate an arbitrary function object over the a set of particles stored
 * in the device. Results are returned to the __host__ via a given mc_host_vector. Datasets with up to nine particles
 * can be handled.
 */
template<typename CUSTOMFUNC, typename RESULT>
void Evaluate(const CUSTOMFUNC funcObj, ParticlesSet_d &pset,
		mc_host_vector<RESULT> &eval)
{

	if (pset.size() > 10 || pset.size() < 2)
	{
		cout
				<< "Can not Calculate(Eval) more than a nine-particle invariant mass."
				<< endl;
		return;
	}

	mc_device_vector<RESULT> dev_out(eval.begin(), eval.end());

	if (pset.size() == 2)
	{
		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];

		thrust::transform(
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end())),
				dev_out.begin(), Calculate<CUSTOMFUNC, RESULT>(funcObj));

	}

	if (pset.size() == 3)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];

		thrust::transform(
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin(),
								dev_v2->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end(),
								dev_v2->end())), dev_out.begin(),
				Calculate<CUSTOMFUNC, RESULT>(funcObj));

	}

	if (pset.size() == 4)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];
		mc_device_vector<Vector4R> *dev_v3 = pset[3];

		thrust::transform(
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin(),
								dev_v2->begin(), dev_v3->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end(),
								dev_v2->end(), dev_v3->end())), dev_out.begin(),
				Calculate<CUSTOMFUNC, RESULT>(funcObj));

	}

	if (pset.size() == 5)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];
		mc_device_vector<Vector4R> *dev_v3 = pset[3];
		mc_device_vector<Vector4R> *dev_v4 = pset[4];

		thrust::transform(
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin(),
								dev_v2->begin(), dev_v3->begin(),
								dev_v4->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end(),
								dev_v2->end(), dev_v3->end(), dev_v4->end())),
				dev_out.begin(), Calculate<CUSTOMFUNC, RESULT>(funcObj));

	}

	if (pset.size() == 6)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];
		mc_device_vector<Vector4R> *dev_v3 = pset[3];
		mc_device_vector<Vector4R> *dev_v4 = pset[4];
		mc_device_vector<Vector4R> *dev_v5 = pset[5];

		thrust::transform(
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin(),
								dev_v2->begin(), dev_v3->begin(),
								dev_v4->begin(), dev_v5->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end(),
								dev_v2->end(), dev_v3->end(), dev_v4->end(),
								dev_v5->end())), dev_out.begin(),
				Calculate<CUSTOMFUNC, RESULT>(funcObj));

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

		thrust::transform(
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin(),
								dev_v2->begin(), dev_v3->begin(),
								dev_v4->begin(), dev_v5->begin(),
								dev_v6->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end(),
								dev_v2->end(), dev_v3->end(), dev_v4->end(),
								dev_v5->end(), dev_v6->end())), dev_out.begin(),
				Calculate<CUSTOMFUNC, RESULT>(funcObj));

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

		thrust::transform(
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin(),
								dev_v2->begin(), dev_v3->begin(),
								dev_v4->begin(), dev_v5->begin(),
								dev_v6->begin(), dev_v7->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end(),
								dev_v2->end(), dev_v3->end(), dev_v4->end(),
								dev_v5->end(), dev_v6->end(), dev_v7->end())),
				dev_out.begin(), Calculate<CUSTOMFUNC, RESULT>(funcObj));

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

		thrust::transform(
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin(),
								dev_v2->begin(), dev_v3->begin(),
								dev_v4->begin(), dev_v5->begin(),
								dev_v6->begin(), dev_v7->begin(),
								dev_v8->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end(),
								dev_v2->end(), dev_v3->end(), dev_v4->end(),
								dev_v5->end(), dev_v6->end(), dev_v7->end(),
								dev_v8->end())), dev_out.begin(),
				Calculate<CUSTOMFUNC, RESULT>(funcObj));

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
		mc_device_vector<Vector4R> *dev_v9 = pset[9];

		thrust::transform(
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin(),
								dev_v2->begin(), dev_v3->begin(),
								dev_v4->begin(), dev_v5->begin(),
								dev_v6->begin(), dev_v7->begin(),
								dev_v8->begin(), dev_v9->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end(),
								dev_v2->end(), dev_v3->end(), dev_v4->end(),
								dev_v5->end(), dev_v6->end(), dev_v7->end(),
								dev_v8->end(), dev_v9->end())), dev_out.begin(),
				Calculate<CUSTOMFUNC, RESULT>(funcObj));

	}

	thrust::copy(dev_out.begin(), dev_out.end(), eval.begin());

	return;
}

/** Template functor for evaluate an arbitrary function object.
 * Template functor for evaluate an arbitrary function object over the a set of particles stored
 * in the device. Results are returned to the __device__ via a given mc_device_vector. Datasets with up to nine particles
 * can be handled.
 */
template<typename CUSTOMFUNC, typename RESULT>
void Evaluate(const CUSTOMFUNC funcObj, ParticlesSet_d &pset,
		mc_device_vector<RESULT> dev_out)
{

	if (pset.size() > 10 || pset.size() < 2)
	{
		cout
				<< "Can not Calculate(Eval) more than a nine-particle invariant mass."
				<< endl;
		return;
	}

	if (pset.size() == 2)
	{
		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];

		thrust::transform(
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end())),
				dev_out.begin(), Calculate<CUSTOMFUNC, RESULT>(funcObj));

	}

	if (pset.size() == 3)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];

		thrust::transform(
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin(),
								dev_v2->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end(),
								dev_v2->end())), dev_out.begin(),
				Calculate<CUSTOMFUNC, RESULT>(funcObj));

	}

	if (pset.size() == 4)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];
		mc_device_vector<Vector4R> *dev_v3 = pset[3];

		thrust::transform(
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin(),
								dev_v2->begin(), dev_v3->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end(),
								dev_v2->end(), dev_v3->end())), dev_out.begin(),
				Calculate<CUSTOMFUNC, RESULT>(funcObj));

	}

	if (pset.size() == 5)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];
		mc_device_vector<Vector4R> *dev_v3 = pset[3];
		mc_device_vector<Vector4R> *dev_v4 = pset[4];

		thrust::transform(
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin(),
								dev_v2->begin(), dev_v3->begin(),
								dev_v4->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end(),
								dev_v2->end(), dev_v3->end(), dev_v4->end())),
				dev_out.begin(), Calculate<CUSTOMFUNC, RESULT>(funcObj));

	}

	if (pset.size() == 6)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];
		mc_device_vector<Vector4R> *dev_v3 = pset[3];
		mc_device_vector<Vector4R> *dev_v4 = pset[4];
		mc_device_vector<Vector4R> *dev_v5 = pset[5];

		thrust::transform(
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin(),
								dev_v2->begin(), dev_v3->begin(),
								dev_v4->begin(), dev_v5->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end(),
								dev_v2->end(), dev_v3->end(), dev_v4->end(),
								dev_v5->end())), dev_out.begin(),
				Calculate<CUSTOMFUNC, RESULT>(funcObj));

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

		thrust::transform(
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin(),
								dev_v2->begin(), dev_v3->begin(),
								dev_v4->begin(), dev_v5->begin(),
								dev_v6->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end(),
								dev_v2->end(), dev_v3->end(), dev_v4->end(),
								dev_v5->end(), dev_v6->end())), dev_out.begin(),
				Calculate<CUSTOMFUNC, RESULT>(funcObj));

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

		thrust::transform(
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin(),
								dev_v2->begin(), dev_v3->begin(),
								dev_v4->begin(), dev_v5->begin(),
								dev_v6->begin(), dev_v7->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end(),
								dev_v2->end(), dev_v3->end(), dev_v4->end(),
								dev_v5->end(), dev_v6->end(), dev_v7->end())),
				dev_out.begin(), Calculate<CUSTOMFUNC, RESULT>(funcObj));

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

		thrust::transform(
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin(),
								dev_v2->begin(), dev_v3->begin(),
								dev_v4->begin(), dev_v5->begin(),
								dev_v6->begin(), dev_v7->begin(),
								dev_v8->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end(),
								dev_v2->end(), dev_v3->end(), dev_v4->end(),
								dev_v5->end(), dev_v6->end(), dev_v7->end(),
								dev_v8->end())), dev_out.begin(),
				Calculate<CUSTOMFUNC, RESULT>(funcObj));

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
		mc_device_vector<Vector4R> *dev_v9 = pset[9];

		thrust::transform(
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin(),
								dev_v2->begin(), dev_v3->begin(),
								dev_v4->begin(), dev_v5->begin(),
								dev_v6->begin(), dev_v7->begin(),
								dev_v8->begin(), dev_v9->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end(),
								dev_v2->end(), dev_v3->end(), dev_v4->end(),
								dev_v5->end(), dev_v6->end(), dev_v7->end(),
								dev_v8->end(), dev_v9->end())), dev_out.begin(),
				Calculate<CUSTOMFUNC, RESULT>(funcObj));

	}

	return;
}

/** Template functor for evaluate an arbitrary function object.
 * Template functor for evaluate an arbitrary function object over the a set of particles stored
 * in the device. No results are returned. This function is used to modify the particles in the input dataset.
 * Datasets with up to nine particles  can be handled.
 */
template<typename CUSTOMFUNC>
void Evaluate(const CUSTOMFUNC funcObj, ParticlesSet_d &pset)
{

	if (pset.size() > 10 || pset.size() < 2)
	{
		cout
				<< "Can not Calculate(Eval) more than a nine-particle invariant mass."
				<< endl;
		return;
	}

	if (pset.size() == 2)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];

		thrust::for_each(
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end())),
				Calculate2<CUSTOMFUNC>(funcObj));

	}

	if (pset.size() == 3)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];

		thrust::for_each(
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin(),
								dev_v2->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end(),
								dev_v2->end())),
				Calculate2<CUSTOMFUNC>(funcObj));

	}

	if (pset.size() == 4)
	{

		mc_device_vector<Vector4R> *dev_v0 = pset[0];
		mc_device_vector<Vector4R> *dev_v1 = pset[1];
		mc_device_vector<Vector4R> *dev_v2 = pset[2];
		mc_device_vector<Vector4R> *dev_v3 = pset[3];

		thrust::for_each(
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin(),
								dev_v2->begin(), dev_v3->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end(),
								dev_v2->end(), dev_v3->end())),
				Calculate2<CUSTOMFUNC>(funcObj));

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
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin(),
								dev_v2->begin(), dev_v3->begin(),
								dev_v4->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end(),
								dev_v2->end(), dev_v3->end(), dev_v4->end())),
				Calculate2<CUSTOMFUNC>(funcObj));

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
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin(),
								dev_v2->begin(), dev_v3->begin(),
								dev_v4->begin(), dev_v5->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end(),
								dev_v2->end(), dev_v3->end(), dev_v4->end(),
								dev_v5->end())),
				Calculate2<CUSTOMFUNC>(funcObj));

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
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin(),
								dev_v2->begin(), dev_v3->begin(),
								dev_v4->begin(), dev_v5->begin(),
								dev_v6->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end(),
								dev_v2->end(), dev_v3->end(), dev_v4->end(),
								dev_v5->end(), dev_v6->end())),
				Calculate2<CUSTOMFUNC>(funcObj));

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
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin(),
								dev_v2->begin(), dev_v3->begin(),
								dev_v4->begin(), dev_v5->begin(),
								dev_v6->begin(), dev_v7->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end(),
								dev_v2->end(), dev_v3->end(), dev_v4->end(),
								dev_v5->end(), dev_v6->end(), dev_v7->end())),
				Calculate2<CUSTOMFUNC>(funcObj));

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
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin(),
								dev_v2->begin(), dev_v3->begin(),
								dev_v4->begin(), dev_v5->begin(),
								dev_v6->begin(), dev_v7->begin(),
								dev_v8->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end(),
								dev_v2->end(), dev_v3->end(), dev_v4->end(),
								dev_v5->end(), dev_v6->end(), dev_v7->end(),
								dev_v8->end())),
				Calculate2<CUSTOMFUNC>(funcObj));

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
		mc_device_vector<Vector4R> *dev_v9 = pset[9];

		thrust::for_each(
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->begin(), dev_v1->begin(),
								dev_v2->begin(), dev_v3->begin(),
								dev_v4->begin(), dev_v5->begin(),
								dev_v6->begin(), dev_v7->begin(),
								dev_v8->begin(), dev_v9->begin())),
				thrust::make_zip_iterator(
						thrust::make_tuple(dev_v0->end(), dev_v1->end(),
								dev_v2->end(), dev_v3->end(), dev_v4->end(),
								dev_v5->end(), dev_v6->end(), dev_v7->end(),
								dev_v8->end(), dev_v9->end())),
				Calculate2<CUSTOMFUNC>(funcObj));

	}

	return;
}

}

#endif /* EVALUATE_H_ */
