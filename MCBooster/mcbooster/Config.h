/*
 * Config.h
 *
 *  Created on: Feb 24, 2016
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

#ifndef CONFIG_H_
#define CONFIG_H_

#include <mcbooster/MCBooster.h>
#include <iostream>

#define CUDA 1
#define OMP 2
#define TBB 3

#if (__cplusplus < 201103L)
#error "This library needs a C++11 compliant compiler"
#endif


#ifndef MCBOOSTER_BACKEND
  #define MCBOOSTER_BACKEND CUDA
#endif


#if MCBOOSTER_BACKEND!=CUDA && MCBOOSTER_BACKEND!=OMP && MCBOOSTER_BACKEND!=TBB

#error "MCBooster: Backend not supported. MCBOOSTER_BACKEND = CUDA, OMP or TBB "

#endif



#if MCBOOSTER_BACKEND==CUDA
	#define CUDA_API_PER_THREAD_DEFAULT_STREAM
	// #define THRUST_DEVICE_SYSTEM CUDA//THRUST_DEVICE_SYSTEM_CUDA
	// #define THRUST_HOST_SYSTEM OMP
	#include <cuda.h>
	#include <cuda_runtime.h>
	#include <cuda_runtime_api.h>
#elif MCBOOSTER_BACKEND==OMP
//    #define THRUST_DEVICE_SYSTEM OMP
//    #define THRUST_HOST_SYSTEM OMP
#elif MCBOOSTER_BACKEND==TBB
//    #define THRUST_DEVICE_SYSTEM TBB
//    #define THRUST_HOST_SYSTEM TBB
#endif


#include <thrust/detail/config/host_device.h>

#endif /* CONFIG_H_ */
