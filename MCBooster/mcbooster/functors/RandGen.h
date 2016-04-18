/*
 * RandGen.cuh
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


#ifndef RANDGEN_H_
#define RANDGEN_H_

#include <mcbooster/Config.h>
#include <thrust/random.h>
#include <mcbooster/GTypes.h>
#include <mcbooster/GContainers.h>


namespace MCBooster
{
/**\struct RandGen
 * Fill a given vector with random numbers between 0 and 1.
 */
struct RandGen
{
	GInt_t fNDaughters; ///< Number of daughter particles
	GReal_t *fRndNumbers;///< Pointer to the array of random numbers
	/**
	 * RandGen ctor. Takes the number of daughter particles and the address of the array
	 * of to be filled with random numbers
	 */
	RandGen(const GInt_t _ndaughters, GReal_t *_rnd) :
		fRndNumbers(_rnd), fNDaughters(_ndaughters)
	{
	}

	/**
	 * hash function. Generate hashs to be used in random number generation initialization
	 */
	__host__      __device__ GUInt_t hash(GUInt_t a)
	{
		a = (a + 0x7ed55d16) + (a << 12);
		a = (a ^ 0xc761c23c) ^ (a >> 19);
		a = (a + 0x165667b1) + (a << 5);
		a = (a + 0xd3a2646c) ^ (a << 9);
		a = (a + 0xfd7046c5) + (a << 3);
		a = (a ^ 0xb55a4f09) ^ (a >> 16);
		return a;
	}

	/**
	 * operator(). Calculate and set random numbers. It takes the index of the event.
	 */
	__host__ __device__ void operator ()(GLong_t idx)
	{
		GUInt_t seed = hash(idx);
     	thrust::random::default_random_engine randEng(seed);
		thrust::uniform_real_distribution<GReal_t> uniDist(0.0, 1.0);


	    fRndNumbers[idx] = uniDist(randEng);

	}

};

struct RandGen2
{
	/**
	 * RandGen2 ctor. Takes the number of daughter particles and the address of the array
	 * of to be filled with random numbers
	 */


	/**
	 * operator(). Calculate and set random numbers. It takes the index of the event.
	 */
	__host__ __device__ GReal_t operator ()(GInt_t idx)
	{

		thrust::random::default_random_engine randEng;
		randEng.discard(idx);
		thrust::uniform_real_distribution<GReal_t> uniDist(0.0, 1.0);

		return uniDist(randEng);


	}

};

}

#endif /* RANDGEN_H_ */
