/*
 * FlagAcceptReject.h
 *
 * Copyright 2016 Antonio Augusto Alves Junior
 *  
 * Created on : 29/03/2016
 *      Author: Antonio Augusto Alves Junior
 */
 
/*
    This file is part of MCBooster.

    MCBooster is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MCBooster is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MCBooster.  If not, see <http://www.gnu.org/licenses/>.
*/

/**\file FlagAcceptReject.h
 * Implements FlagAcceptReject.
 */
#ifndef FLAGACCEPTEDREJECTED_H_
#define FLAGACCEPTEDREJECTED_H_



#include <mcbooster/Config.h>
#include <thrust/random.h>
#include <mcbooster/GTypes.h>

namespace MCBooster
{
/**\struct FlagAcceptReject
 * Flags generated events as accepted (1) or rejected (0).
 */
struct FlagAcceptReject
{

	GReal_t wmax; ///< maximum weight
	/**
	 * FlagAcceptReject constructor. It is initialized with the value of the maximum weight
	 * with which event weights will be compared.
	 */
	FlagAcceptReject(const GReal_t _wmax) :
		wmax(_wmax)
	{	}

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
	 * operator(). Takes the events index and weight and so flag it as accepted and rejected
	 *
	 */
	__host__ __device__ GBool_t operator ()(GLong_t idx, GReal_t weight)
	{
		GUInt_t seed = hash(idx+68464654684);
		thrust::default_random_engine randEng(seed);
		thrust::uniform_real_distribution<GReal_t> uniDist(0.0, wmax);


		GBool_t flag = (uniDist(randEng) < weight) ? 1 : 0;
		return flag;

	}

};
}



#endif /* FLAGACCEPTEDREJECTED_H_ */
