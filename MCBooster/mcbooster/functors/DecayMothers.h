/*
 * DecayMothers.cuh
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

#ifndef DECAYMOTHERS_H_
#define DECAYMOTHERS_H_


#include <mcbooster/Config.h>
#include <mcbooster/GTypes.h>
#include <mcbooster/GContainers.h>
#include <mcbooster/Vector3R.h>
#include <mcbooster/Vector4R.h>
#include <thrust/tuple.h>
#include <thrust/iterator/zip_iterator.h>

using namespace std;

namespace MCBooster
{

struct DecayMothers
{
	const GInt_t fSeed;
	const GInt_t fNDaughters;
	const GReal_t* __restrict__ fMasses;

	//constructor
	DecayMothers(const mc_device_vector<GReal_t>& _masses,
			const GInt_t _ndaughters, const GInt_t _seed ):
			fMasses(thrust::raw_pointer_cast(_masses.data())),
			fNDaughters(_ndaughters),
			fSeed(_seed)
	{
	}

	__host__      __device__ GReal_t pdk(const GReal_t a, const GReal_t b,
			const GReal_t c) const
	{
		//the PDK function
		GReal_t x = (a - b - c) * (a + b + c) * (a - b + c) * (a + b - c);
		x = sqrt(x) / (2 * a);
		return x;
	}

	__host__ __device__ void bbsort(GReal_t *array, GInt_t n)
	{
		// Improved bubble sort
		for (GInt_t c = 0; c < n; c++)
		{
			GInt_t nswap = 0;

			for (GInt_t d = 0; d < n - c - 1; d++)
			{
				if (array[d] > array[d + 1]) /* For decreasing order use < */
				{
					GReal_t swap = array[d];
					array[d] = array[d + 1];
					array[d + 1] = swap;
					nswap++;
				}
			}
			if (nswap == 0)
				break;
		}

	}

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

	__host__      __device__ GReal_t process(const GInt_t evt, Vector4R** particles)
	{


		thrust::random::default_random_engine randEng( hash(evt)*fSeed);
		thrust::uniform_real_distribution<GReal_t> uniDist(0.0, 1.0);

		GReal_t fTeCmTm = 0.0, fWtMax = 0.0;

		fTeCmTm = particles[0]->mass(); // total energy in C.M. minus the sum of the masses

		#pragma unroll 9
		for (size_t n = 0; n < fNDaughters; n++)
		{
			fTeCmTm -= fMasses[n];
		}

		GReal_t emmax = fTeCmTm + fMasses[0];
		GReal_t emmin = 0.0;
		GReal_t wtmax = 1.0;

		#pragma unroll 9
		for (size_t n = 1; n < fNDaughters; n++)
		{
			emmin += fMasses[n - 1];
			emmax += fMasses[n];
			wtmax *= pdk(emmax, emmin, fMasses[n]);
		}
		fWtMax = 1.0 / wtmax;
		//
		//---->  get the betas of the decaying particle
		//
		GReal_t fBeta[3];
		fBeta[0]=0, fBeta[1]=0, fBeta[2] = 0.0;

		GReal_t _beta = particles[0]->d3mag() / particles[0]->get(0);

		if (_beta)
		{
			GReal_t w = _beta / particles[0]->d3mag();
			fBeta[0] = particles[0]->get(0) * w;
			fBeta[1] = particles[0]->get(1) * w;
			fBeta[2] = particles[0]->get(2) * w;
		}

		GReal_t rno[kMAXP];
		rno[0] = 0.0;

		if (fNDaughters > 2)
		{
			#pragma unroll 9
			for (size_t n = 1; n < fNDaughters - 1; n++)
				rno[n] = uniDist(randEng) ;
			bbsort(&rno[1], fNDaughters - 2);

		}
		rno[fNDaughters - 1] = 1;

		GReal_t invMas[kMAXP], sum = 0.0;

		#pragma unroll 9
		for (size_t n = 0; n < fNDaughters; n++)
		{
			sum += fMasses[n];
			invMas[n] = rno[n] * fTeCmTm + sum;
		}

		//
		//-----> compute the weight of the current event
		//

		GReal_t wt  = 1.0 / wtmax;

		GReal_t pd[kMAXP];

		#pragma unroll 9
		for (size_t n = 0; n < fNDaughters - 1; n++)
		{
			pd[n] = pdk(invMas[n + 1], invMas[n], fMasses[n + 1]);
			wt *= pd[n];
		}

		//
		//-----> complete specification of event (Raubold-Lynch method)
		//

		particles[1]->set(sqrt(pd[0] * pd[0] + fMasses[0] * fMasses[0]), 0.0,
				pd[0], 0.0);

		#pragma unroll 9
		for (size_t i = 1; i < fNDaughters; i++)
		{

			particles[i + 1]->set(
					sqrt(pd[i - 1] * pd[i - 1] + fMasses[i] * fMasses[i]), 0.0,
					-pd[i - 1], 0.0);

			GReal_t cZ = 2	* uniDist(randEng) -1 ;
			GReal_t sZ = sqrt(1 - cZ * cZ);
			GReal_t angY = 2.0 * PI	* uniDist(randEng);
			GReal_t cY = cos(angY);
			GReal_t sY = sin(angY);
			for (size_t j = 0; j <= i; j++)
			{

				GReal_t x = particles[j + 1]->get(1);
				GReal_t y = particles[j + 1]->get(2);
				particles[j + 1]->set(1, cZ * x - sZ * y);
				particles[j + 1]->set(2, sZ * x + cZ * y); // rotation around Z

				x = particles[j + 1]->get(1);
				GReal_t z = particles[j + 1]->get(3);
				particles[j + 1]->set(1, cY * x - sY * z);
				particles[j + 1]->set(3, sY * x + cY * z); // rotation around Y
			}

			if (i == (fNDaughters - 1))
				break;

			GReal_t beta = pd[i] / sqrt(pd[i] * pd[i] + invMas[i] * invMas[i]);
			for (size_t j = 0; j <= i; j++)
			{

				particles[j + 1]->applyBoostTo(0, beta, 0);
			}

			//i++;
		}

		//
		//---> final boost of all particles to the mother's frame
		//
		#pragma unroll 9
		for (size_t n = 0; n < fNDaughters; n++)
		{

			particles[n + 1]->applyBoostTo(*particles[0]);

		}

		//
		//---> return the weight of event
		//

		return wt;

	}

	__host__      __device__ GReal_t operator()(const GInt_t evt, GT2 &particles)
	{
		//do nothing, will never be called
		return 0.0;

	}

	__host__      __device__ GReal_t operator()(const GInt_t evt, GT3& particles)
	{
		Vector4R* _Particles[3];

		_Particles[0] = &thrust::get<0>(particles);
		_Particles[1] = &thrust::get<1>(particles);
		_Particles[2] = &thrust::get<2>(particles);

		return process(evt, _Particles);

	}

	__host__      __device__ GReal_t operator()(const GInt_t evt, GT4& particles)
	{

		Vector4R* _Particles[4];

		_Particles[0] = &thrust::get<0>(particles);
		_Particles[1] = &thrust::get<1>(particles);
		_Particles[2] = &thrust::get<2>(particles);
		_Particles[3] = &thrust::get<3>(particles);

		return process(evt, _Particles);

	}

	__host__      __device__ GReal_t operator()(const GInt_t evt, GT5& particles)
	{
		Vector4R* _Particles[5];

		_Particles[0] = &thrust::get<0>(particles);
		_Particles[1] = &thrust::get<1>(particles);
		_Particles[2] = &thrust::get<2>(particles);
		_Particles[3] = &thrust::get<3>(particles);
		_Particles[4] = &thrust::get<4>(particles);

		return process(evt, _Particles);

	}

	__host__      __device__ GReal_t operator()(const GInt_t evt, GT6& particles)
	{
		Vector4R* _Particles[6];

		_Particles[0] = &thrust::get<0>(particles);
		_Particles[1] = &thrust::get<1>(particles);
		_Particles[2] = &thrust::get<2>(particles);
		_Particles[3] = &thrust::get<3>(particles);
		_Particles[4] = &thrust::get<4>(particles);
		_Particles[5] = &thrust::get<5>(particles);

		return process(evt, _Particles);

	}

	__host__      __device__ GReal_t operator()(const GInt_t evt, GT7& particles)
	{
		Vector4R* _Particles[7];

		_Particles[0] = &thrust::get<0>(particles);
		_Particles[1] = &thrust::get<1>(particles);
		_Particles[2] = &thrust::get<2>(particles);
		_Particles[3] = &thrust::get<3>(particles);
		_Particles[4] = &thrust::get<4>(particles);
		_Particles[5] = &thrust::get<5>(particles);
		_Particles[6] = &thrust::get<6>(particles);

		return process(evt, _Particles);

	}

	__host__      __device__ GReal_t operator()(const GInt_t evt, GT8& particles)
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

		return process(evt, _Particles);

	}

	__host__      __device__ GReal_t operator()(const GInt_t evt, GT9& particles)
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

		return process(evt, _Particles);

	}

	__host__      __device__ GReal_t operator()(const GInt_t evt, GT10& particles)
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

		return process(evt, _Particles);

	}

};

}

#endif /* DECAYMOTHERS_H_ */
