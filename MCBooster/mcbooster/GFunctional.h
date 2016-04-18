/*
 * GFunctional.h
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

/** \file GFunctional.h
 * Implements the template functors IFunction and IFunctionArray
 */
#ifndef GFUNCTIONAL_H_
#define GFUNCTIONAL_H_


#include <mcbooster/Config.h>
#include <mcbooster/GTypes.h>
#include <mcbooster/Vector3R.h>
#include <mcbooster/Vector4R.h>

namespace MCBooster
{
/** \struct  IFunction
 *  IFunction is the base class for arbitrary functions return any type suported by the framwork.
 */
template<typename RESULT>
struct IFunction
{

	__host__   __device__   virtual RESULT operator()(const GInt_t n,
			Vector4R** particles)=0;

};

/** \struct  IFunction
 *  IFunctionArray is the base class for arbitrary functions used to evaluate at once an array of variables.
 */
struct IFunctionArray
{
	GInt_t dim;
	IFunctionArray() :
			dim(0)
	{	}

	__host__ __device__ virtual void operator()(const GInt_t np,
			Vector4R** particles, GReal_t* variables)=0;

};

}

#endif /* GFUNCTIONAL_H_ */
