/*
 * Calculate.cuh
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

#ifndef CALCULATE_H_
#define CALCULATE_H_

#include <mcbooster/Config.h>
#include <mcbooster/Vector3R.h>
#include <mcbooster/Vector4R.h>
#include <mcbooster/GContainers.h>
#include <mcbooster/GTypes.h>

namespace MCBooster
{

template<typename FUNCTION, typename RESULT>
struct Calculate
{

	FUNCTION Function;

	Calculate()
	{

		Function = FUNCTION();

	}


	Calculate(const FUNCTION& _Function) :
			Function(_Function)
	{

	}


	__host__      __device__ RESULT operator()(GT2 &particles)
	{

		Vector4R* _Particles[2];

		_Particles[0] = &thrust::get<0>(particles);
		_Particles[1] = &thrust::get<1>(particles);

		return Function(2, _Particles);
	}

	__host__      __device__ RESULT operator()(GT3 &particles)
	{

		Vector4R* _Particles[3];

		_Particles[0] = &thrust::get<0>(particles);
		_Particles[1] = &thrust::get<1>(particles);
		_Particles[2] = &thrust::get<2>(particles);

		return Function(3, _Particles);
	}

	__host__         __device__ RESULT operator()(GT4 &particles)
	{

		Vector4R* _Particles[4];

		_Particles[0] = &thrust::get<0>(particles);
		_Particles[1] = &thrust::get<1>(particles);
		_Particles[2] = &thrust::get<2>(particles);
		_Particles[3] = &thrust::get<3>(particles);

		return Function(4, _Particles);

	}

	__host__         __device__ RESULT operator()(GT5 &particles)
	{

		Vector4R* _Particles[5];

		_Particles[0] = &thrust::get<0>(particles);
		_Particles[1] = &thrust::get<1>(particles);
		_Particles[2] = &thrust::get<2>(particles);
		_Particles[3] = &thrust::get<3>(particles);
		_Particles[4] = &thrust::get<4>(particles);

		return Function(5, _Particles);
	}

	__host__         __device__ RESULT operator()(GT6 &particles)
	{

		Vector4R* _Particles[6];

		_Particles[0] = &thrust::get<0>(particles);
		_Particles[1] = &thrust::get<1>(particles);
		_Particles[2] = &thrust::get<2>(particles);
		_Particles[3] = &thrust::get<3>(particles);
		_Particles[4] = &thrust::get<4>(particles);
		_Particles[5] = &thrust::get<5>(particles);

		return Function(6, _Particles);
	}

	__host__         __device__ RESULT operator()(GT7 &particles)
	{

		Vector4R* _Particles[7];

		_Particles[0] = &thrust::get<0>(particles);
		_Particles[1] = &thrust::get<1>(particles);
		_Particles[2] = &thrust::get<2>(particles);
		_Particles[3] = &thrust::get<3>(particles);
		_Particles[4] = &thrust::get<4>(particles);
		_Particles[5] = &thrust::get<5>(particles);
		_Particles[6] = &thrust::get<6>(particles);

		return Function(7, _Particles);
	}

	__host__         __device__ RESULT operator()(GT8 &particles)
	{

		Vector4R* _Particles[8];

		_Particles[0] = &thrust::get<0>(particles);
		_Particles[1] = &thrust::get<1>(particles);
		_Particles[2] = &thrust::get<2>(particles);
		_Particles[3] = &thrust::get<3>(particles);
		_Particles[4] = &thrust::get<4>(particles);
		_Particles[5] = &thrust::get<5>(particles);
		_Particles[6] = &thrust::get<6>(particles);
		_Particles[7] = &thrust::get<7>(particles);

		return Function(8, _Particles);

	}

	__host__         __device__ RESULT operator()(GT9 &particles)
	{

		Vector4R* _Particles[9];

		_Particles[0] = &thrust::get<0>(particles);
		_Particles[1] = &thrust::get<1>(particles);
		_Particles[2] = &thrust::get<2>(particles);
		_Particles[3] = &thrust::get<3>(particles);
		_Particles[4] = &thrust::get<4>(particles);
		_Particles[5] = &thrust::get<5>(particles);
		_Particles[6] = &thrust::get<6>(particles);
		_Particles[7] = &thrust::get<7>(particles);
		_Particles[8] = &thrust::get<8>(particles);

		return Function(9, _Particles);
	}

	__host__         __device__ RESULT operator()(GT10 &particles)
	{

		Vector4R* _Particles[10];

		_Particles[0] = &thrust::get<0>(particles);
		_Particles[1] = &thrust::get<1>(particles);
		_Particles[2] = &thrust::get<2>(particles);
		_Particles[3] = &thrust::get<3>(particles);
		_Particles[4] = &thrust::get<4>(particles);
		_Particles[5] = &thrust::get<5>(particles);
		_Particles[6] = &thrust::get<6>(particles);
		_Particles[7] = &thrust::get<7>(particles);
		_Particles[8] = &thrust::get<8>(particles);
		_Particles[9] = &thrust::get<9>(particles);

		return Function(10, _Particles);
	}

};

template<typename FUNCTION>
struct Calculate2
{

	FUNCTION Function;

	Calculate2()
	{
		Function = FUNCTION();

	}


	Calculate2(const FUNCTION& _Function) :
			Function(_Function)
	{

	}


	__host__ __device__ void operator()(GT2 particles)
	{

		Vector4R* _Particles[2];

		_Particles[0] = &thrust::get<0>(particles);
		_Particles[1] = &thrust::get<1>(particles);

		Function(2, _Particles);
	}

	__host__ __device__ void operator()(GT3 particles)
	{

		Vector4R* _Particles[3];

		_Particles[0] = &thrust::get<0>(particles);
		_Particles[1] = &thrust::get<1>(particles);
		_Particles[2] = &thrust::get<2>(particles);

		Function(3, _Particles);
	}

	__host__ __device__ void operator()(GT4 particles)
	{

		Vector4R* _Particles[4];

		_Particles[0] = &thrust::get<0>(particles);
		_Particles[1] = &thrust::get<1>(particles);
		_Particles[2] = &thrust::get<2>(particles);
		_Particles[3] = &thrust::get<3>(particles);

		Function(4, _Particles);

	}

	__host__ __device__ void operator()(GT5 particles)
	{

		Vector4R* _Particles[5];

		_Particles[0] = &thrust::get<0>(particles);
		_Particles[1] = &thrust::get<1>(particles);
		_Particles[2] = &thrust::get<2>(particles);
		_Particles[3] = &thrust::get<3>(particles);
		_Particles[4] = &thrust::get<4>(particles);

		Function(5, _Particles);
	}

	__host__ __device__ void operator()(GT6 particles)
	{

		Vector4R* _Particles[6];

		_Particles[0] = &thrust::get<0>(particles);
		_Particles[1] = &thrust::get<1>(particles);
		_Particles[2] = &thrust::get<2>(particles);
		_Particles[3] = &thrust::get<3>(particles);
		_Particles[4] = &thrust::get<4>(particles);
		_Particles[5] = &thrust::get<5>(particles);

		Function(6, _Particles);
	}

	__host__ __device__ void operator()(GT7 particles)
	{

		Vector4R* _Particles[7];

		_Particles[0] = &thrust::get<0>(particles);
		_Particles[1] = &thrust::get<1>(particles);
		_Particles[2] = &thrust::get<2>(particles);
		_Particles[3] = &thrust::get<3>(particles);
		_Particles[4] = &thrust::get<4>(particles);
		_Particles[5] = &thrust::get<5>(particles);
		_Particles[6] = &thrust::get<6>(particles);

		Function(7, _Particles);
	}

	__host__ __device__ void operator()(GT8 particles)
	{

		Vector4R* _Particles[8];

		_Particles[0] = &thrust::get<0>(particles);
		_Particles[1] = &thrust::get<1>(particles);
		_Particles[2] = &thrust::get<2>(particles);
		_Particles[3] = &thrust::get<3>(particles);
		_Particles[4] = &thrust::get<4>(particles);
		_Particles[5] = &thrust::get<5>(particles);
		_Particles[6] = &thrust::get<6>(particles);
		_Particles[7] = &thrust::get<7>(particles);

		Function(8, _Particles);

	}

	__host__ __device__ void operator()(GT9 particles)
	{

		Vector4R* _Particles[9];

		_Particles[0] = &thrust::get<0>(particles);
		_Particles[1] = &thrust::get<1>(particles);
		_Particles[2] = &thrust::get<2>(particles);
		_Particles[3] = &thrust::get<3>(particles);
		_Particles[4] = &thrust::get<4>(particles);
		_Particles[5] = &thrust::get<5>(particles);
		_Particles[6] = &thrust::get<6>(particles);
		_Particles[7] = &thrust::get<7>(particles);
		_Particles[8] = &thrust::get<8>(particles);

		Function(9, _Particles);
	}

	__host__ __device__ void operator()(GT10 particles)
	{

		Vector4R* _Particles[10];

		_Particles[0] = &thrust::get<0>(particles);
		_Particles[1] = &thrust::get<1>(particles);
		_Particles[2] = &thrust::get<2>(particles);
		_Particles[3] = &thrust::get<3>(particles);
		_Particles[4] = &thrust::get<4>(particles);
		_Particles[5] = &thrust::get<5>(particles);
		_Particles[6] = &thrust::get<6>(particles);
		_Particles[7] = &thrust::get<7>(particles);
		_Particles[8] = &thrust::get<8>(particles);
		_Particles[9] = &thrust::get<9>(particles);

		Function(10, _Particles);
	}

};

template<typename FUNCTION>
struct Calculate3
{

	FUNCTION Function;

	Calculate3()
	{
		Function = FUNCTION();

	}


	Calculate3(const FUNCTION& _Function) :
			Function(_Function)
	{

	}


	__host__ __device__ void operator()(GTR3 tuples)
	{

		GReal_t* _real;
		Vector4R* _Particles[2];

		_real = &thrust::get<0>(tuples);
		_Particles[0] = &thrust::get<1>(tuples);
		_Particles[1] = &thrust::get<2>(tuples);

		Function(2, _Particles, _real);
	}

	__host__ __device__ void operator()(GTR4 tuples)
	{

		GReal_t* _real;

		Vector4R* _Particles[3];

		_real = &thrust::get<0>(tuples);
		_Particles[0] = &thrust::get<1>(tuples);
		_Particles[1] = &thrust::get<2>(tuples);
		_Particles[2] = &thrust::get<3>(tuples);

		Function(3, _Particles, _real);

	}

	__host__ __device__ void operator()(GTR5 tuples)
	{

		GReal_t* _real;

		Vector4R* _Particles[4];

		_real = &thrust::get<0>(tuples);
		_Particles[0] = &thrust::get<1>(tuples);
		_Particles[1] = &thrust::get<2>(tuples);
		_Particles[2] = &thrust::get<3>(tuples);
		_Particles[3] = &thrust::get<4>(tuples);

		Function(4, _Particles, _real);
	}

	__host__ __device__ void operator()(GTR6 tuples)
	{

		GReal_t* _real;

		Vector4R* _Particles[5];

		_real = &thrust::get<0>(tuples);
		_Particles[0] = &thrust::get<1>(tuples);
		_Particles[1] = &thrust::get<2>(tuples);
		_Particles[2] = &thrust::get<3>(tuples);
		_Particles[3] = &thrust::get<4>(tuples);
		_Particles[4] = &thrust::get<5>(tuples);

		Function(5, _Particles, _real);
	}

	__host__ __device__ void operator()(GTR7 tuples)
	{

		GReal_t* _real;

		Vector4R* _Particles[6];

		_real = &thrust::get<0>(tuples);
		_Particles[0] = &thrust::get<1>(tuples);
		_Particles[1] = &thrust::get<2>(tuples);
		_Particles[2] = &thrust::get<3>(tuples);
		_Particles[3] = &thrust::get<4>(tuples);
		_Particles[4] = &thrust::get<5>(tuples);
		_Particles[5] = &thrust::get<6>(tuples);

		Function(6, _Particles, _real);
	}

	__host__ __device__ void operator()(GTR8 tuples)
	{

		GReal_t* _real;

		Vector4R* _Particles[7];

		_real = &thrust::get<0>(tuples);
		_Particles[0] = &thrust::get<1>(tuples);
		_Particles[1] = &thrust::get<2>(tuples);
		_Particles[2] = &thrust::get<3>(tuples);
		_Particles[3] = &thrust::get<4>(tuples);
		_Particles[4] = &thrust::get<5>(tuples);
		_Particles[5] = &thrust::get<6>(tuples);
		_Particles[6] = &thrust::get<7>(tuples);

		Function(7, _Particles, _real);

	}

	__host__ __device__ void operator()(GTR9 tuples)
	{

		GReal_t* _real;

		Vector4R* _Particles[8];

		_real = &thrust::get<0>(tuples);
		_Particles[0] = &thrust::get<1>(tuples);
		_Particles[1] = &thrust::get<2>(tuples);
		_Particles[2] = &thrust::get<3>(tuples);
		_Particles[3] = &thrust::get<4>(tuples);
		_Particles[4] = &thrust::get<5>(tuples);
		_Particles[5] = &thrust::get<6>(tuples);
		_Particles[6] = &thrust::get<7>(tuples);
		_Particles[7] = &thrust::get<8>(tuples);

		Function(8, _Particles, _real);
	}

	__host__ __device__ void operator()(GTR10 tuples)
	{

		GReal_t* _real;

		Vector4R* _Particles[9];

		_real = &thrust::get<0>(tuples);
		_Particles[0] = &thrust::get<1>(tuples);
		_Particles[1] = &thrust::get<2>(tuples);
		_Particles[2] = &thrust::get<3>(tuples);
		_Particles[3] = &thrust::get<4>(tuples);
		_Particles[4] = &thrust::get<5>(tuples);
		_Particles[5] = &thrust::get<6>(tuples);
		_Particles[6] = &thrust::get<7>(tuples);
		_Particles[7] = &thrust::get<8>(tuples);
		_Particles[8] = &thrust::get<9>(tuples);

		Function(9, _Particles, _real);
	}

};
}

#endif /* CALCULATE_H_ */
