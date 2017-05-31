/*-
 * Generate.h
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

/*!\file Generate.h
 * Implements the struct Events and the class PhaseSpace
 */

#ifndef GENERATE_H_
#define GENERATE_H_

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <ostream>
#include <algorithm>
#include <ctime>
#include <cstdio>
//#include <math.h>

#include <mcbooster/Config.h>
#include <mcbooster/Vector3R.h>
#include <mcbooster/Vector4R.h>
#include <mcbooster/GTypes.h>
#include <mcbooster/GContainers.h>
#include <mcbooster/functors/DecayMother.h>
#include <mcbooster/functors/DecayMothers.h>
#include <mcbooster/functors/RandGen.h>
#include <mcbooster/functors/FlagAcceptReject.h>
#include <mcbooster/functors/IsAccepted.h>
#include <mcbooster/strided_iterator.h>

#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/transform.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/sequence.h>
#include <thrust/for_each.h>
#include <thrust/tuple.h>
#include <thrust/extrema.h>
#include <thrust/count.h>
#include <thrust/fill.h>
#include <thrust/sort.h>
#include <thrust/iterator/counting_iterator.h>

#include <thrust/execution_policy.h>

#define TIMER CLOCK_REALTIME

#define CUDA_CHECK_RETURN(value)                                                                                       \
    {                                                                                                                  \
        cudaError_t _m_cudaStat = value;                                                                               \
        if(_m_cudaStat != cudaSuccess) {                                                                               \
            fprintf(stderr, "Error %s at line %d in file %s\n", cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__);  \
            exit(1);                                                                                                   \
        }                                                                                                              \
    }

namespace mcbooster {
/*!
 * Function to calculate time intervals in seconds.
 */
inline timespec time_diff(timespec start, timespec end) {
    timespec temp;
    if((end.tv_nsec - start.tv_nsec) < 0) {
        temp.tv_sec  = end.tv_sec - start.tv_sec - 1;
        temp.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
    } else {
        temp.tv_sec  = end.tv_sec - start.tv_sec;
        temp.tv_nsec = end.tv_nsec - start.tv_nsec;
    }
    return temp;
}

/*! \struct Events
 * Events is a container struct to hold all the information corresponding the generated events.
 * Mother four-vectors are not stored.
 */
struct Events {
    GInt_t fNDaughters;            ///< Number of daughters.
    GLong_t fNEvents;              ///< Number of events.
    GReal_t fMaxWeight;            ///< Maximum weight of the generated events.
    BoolVector_h fAccRejFlags;     ///< Vector of flags. Accepted events are flagged 1 and rejected 0.
    RealVector_h fWeights;         ///< Vector of event weights.
    Particles_h fDaughters[kMAXP]; ///< Array of daughter particle vectors.

    /*!
     * Constructor takes as parameters the number of particles and number of events.
     */
    Events(GInt_t ndaughters, GLong_t nevents)
        : fNDaughters(ndaughters)
        , fNEvents(nevents)
        , fMaxWeight(0) {
        for(GInt_t d = 0; d < fNDaughters; d++) {
            fDaughters[d].resize(fNEvents);
        }

        fWeights.resize(fNEvents);
        fAccRejFlags.resize(fNEvents);
    }
};

/*!\class PhaseSpace
 *Implementation of the event generator.
 *Basic usage:
 *@code
 *#include <mcbooster/GTypes.h>
 *#include <mcbooster/Vector4R.h>
 *#include <mcbooster/Generate.h>
 * ...
 * //setting the mother particle
 * Vector4R B0(5.2795, 0.0, 0.0, 0.0);
 *
 * //setting the masses of the daughter particles
 * vector<GReal_t> masses;
 * masses.push_back(3.096916); // J/psi
 * masses.push_back(0.493677); // K
 * masses.push_back(0.13957018); // pi
 *
 * //generator ctor for 10M events
 * PhaseSpace phsp(B0.mass(), massesB0, 10000000);
 *
 * //run the generator
 * 	  phsp.Generate(B0);
 *
 * //flag the accepted and rejected events
 * phsp.Unweight();
 *
 * //export events to the host (in case it is necessary)
 * Events *GenEvents = new Events(masses.size(), 10000000);
 * phsp.Export(GenEvents);
 * ...
 *@endcode
 */
class PhaseSpace {
  public:
    /**
     * PhaseSpace ctor. Constructor of the phase-space generator takes as input parameters:
     * - _MotherMass: the mass of the mother particle in Gev/c*c
     * - _Masses: STL vector with the mass of the daughter particles.
     * - _NEvents: it is the number of events to be generated.
     */
    PhaseSpace(GReal_t _MotherMass, std::vector<GReal_t> _Masses, GLong_t _NEvents)
        : fNEvents(_NEvents)
        , fNDaughters(_Masses.size())
        , fSeed(0) {
        if(_Masses.size() < 2 || _Masses.size() > 9) {
            std::cout << "The number of daughter particles can not be (< 2 || > 9) or masses and names need to have "
                         "the same size."
                      << std::endl;
            exit(1);
        }

        fMasses.resize(_Masses.size());
        thrust::copy(_Masses.begin(), _Masses.end(), fMasses.begin());

        GReal_t fTeCmTm = 0.0;

        fTeCmTm = _MotherMass; // total energy in C.M. minus the sum of the masses

        for(size_t n = 0; n < fNDaughters; n++) {
            fTeCmTm -= _Masses[n];
        }
        if(fTeCmTm < 0.0) {
            std::cout << "Not enough energy for this decay. Exit." << std::endl;
            exit(1);
        }

        Allocate(fNEvents);

    } // decay

    PhaseSpace(GReal_t _MotherMass, std::vector<GReal_t> _Masses, GLong_t _NEvents, GLong_t _EvtNumOffset)
        : PhaseSpace(_MotherMass, _Masses, _NEvents) {
        fSeed = _EvtNumOffset;
    }
    /**
     * PhaseSpace() dtor. The destructor explicitly frees all the resources owned by the class.
     */
    ~PhaseSpace() {
        fMasses.clear();
        fMasses.shrink_to_fit();
        FreeResources();
    }

    /**
     * Free resources. Free resources owned in host and device side by the PhaseSpace object.
     */
    inline void FreeResources() {
        for(GInt_t i = 0; i < fNDaughters; i++) {
            fDaughters[i].clear();
            fDaughters[i].shrink_to_fit();
        }

        fWeights.clear();
        fWeights.shrink_to_fit();
        // fRandNumbers.clear();
        // fRandNumbers.shrink_to_fit();
        fAccRejFlags.clear();
        fAccRejFlags.shrink_to_fit();
    }

    inline void Generate(Particles_d fMothers);
    inline void Generate(const Vector4R fMother);

    /**
     * Get the daughter with index 'i' in the mass array. It return a device vector of particles by reference.
     * Is responsibility of the user modify or not the particles in the container.
     */
    inline Particles_d &GetDaughters(GInt_t i) { return fDaughters[i]; }
    /**
     * Returns the number of daughter particles.
     */
    inline GInt_t GetNDaughters() const { return fNDaughters; }

    /**
     * Returns the number of events.
     */
    inline GLong_t GetNEvents() const { return fNEvents; }

    inline GLong_t GetNAccepted() const { return fNAccepted; }

    inline BoolVector_d GetAccRejFlags() const { return fAccRejFlags; }

    /**
     * Returns a device vector with the event weights.
     */
    inline const RealVector_d &GetWeights() const { return fWeights; }

    /**
     * Returns the time spent in seconds to process the random numbers and set the phase space four-vectors.
     */
    inline GReal_t GetEvtTime() const { return EVT_Time; }

    /**
     * Returns the time spent in seconds to export the generated events to a Events container.
     */
    inline GReal_t GetExpTime() const { return EXP_Time; }

    /**
     * Returns the time spent in seconds to generate the random numbers necessary for the calculation.
     */
    inline GReal_t GetRndTime() const { return RND_Time; }

    /**
     * Returns the max event weight that can be found in the generated sample. Larger the sample, more accurate
     * is this number.
     */
    inline GReal_t GetMaxWeight() const { return fMaxWeight; }

    inline GInt_t GetSeed() const { return fSeed; }

    inline void SetSeed(GInt_t _seed) { fSeed = _seed; }

    /**
     * Export the events and all related information to host.
     */
    void Export(Events *_Events);
    inline void ExportUnweighted(Events *_Events);
    /**
     * Flag the accepted and rejected events
     */
    inline GULong_t Unweight();

    //
  public:
    /**
     * Allocate resources on the device for event generation.
     */
    inline void Allocate(const GLong_t _nevents) {
        for(GInt_t i = 0; i < fNDaughters; i++)
            fDaughters[i].resize(_nevents);
        fWeights.resize(_nevents);
        fAccRejFlags.resize(_nevents);
        // fRandNumbers.resize((3 * fNDaughters - 2) * fNEvents);
    }

    /**
     * PDK function
     */
    inline GReal_t PDK(const GReal_t a, const GReal_t b, const GReal_t c) const {
        // the PDK function
        GReal_t x = (a - b - c) * (a + b + c) * (a - b + c) * (a + b - c);
        x         = sqrt(x) / (2 * a);
        return x;
    }

  private:
    GLong_t fNEvents;   ///< Number of events.
    GInt_t fNDaughters; ///< Number of daughters.
    GUInt_t fSeed{0};   ///< seed.
    GLong_t fNAccepted;
    GReal_t RND_Time{0.0};   ///< Random number generation time interval seconds.
    GReal_t EVT_Time{0.0};   ///< Event generation time interval in seconds.
    GReal_t EXP_Time{0.0};   ///< Events export time interval in seconds.
    GReal_t fMaxWeight{0.0}; ///< Maximum weight in sample.
    // device
    RealVector_d fMasses;          ///< Device vector of daughter masses.
    RealVector_d fWeights;         ///< Device vector of weights.
    BoolVector_d fAccRejFlags;     ///< Device vector of Accept/reject flags
    Particles_d fDaughters[kMAXP]; ///< Array of device vectors with the daughter four-vectors
};

inline GULong_t PhaseSpace::Unweight() {
    /**
     * Flag the accepted and rejected events
     */

    GULong_t count = 0;
    if(fNDaughters == 2) {
        thrust::fill(fAccRejFlags.begin(), fAccRejFlags.end(), kTrue);
        count = fNEvents;
    } else {
        // create iterators
        thrust::counting_iterator<GLong_t> first(0);
        thrust::counting_iterator<GLong_t> last = first + fNEvents;

        thrust::transform(first, last, fWeights.begin(), fAccRejFlags.begin(), FlagAcceptReject(fMaxWeight, fSeed));

        count = thrust::count(fAccRejFlags.begin(), fAccRejFlags.end(), kTrue);
    }
    fNAccepted = count;
    return count;
}

inline void PhaseSpace::ExportUnweighted(Events *_Events) {
    /**
     * Export the events and all related information to an Events object properly initialized.
     */

    if(!fNAccepted)
        Unweight();

    _Events->fMaxWeight = fMaxWeight;

#if MCBOOSTER_BACKEND != CUDA
#if MCBOOSTER_BACKEND != CPP
#pragma omp parallel num_threads(fNDaughters + 1)
    {
#else
    for(size_t val = 0; val < fNDaughters; val++) {
        auto omp_get_thread_num = [&val]() { return val; };
#endif
        if(omp_get_thread_num() < fNDaughters) {
            thrust::copy_if(fDaughters[omp_get_thread_num()].begin(),
                            fDaughters[omp_get_thread_num()].end(),
                            fAccRejFlags.begin(),
                            _Events->fDaughters[omp_get_thread_num()].begin(),
                            isAccepted());
        }

        if(omp_get_thread_num() == fNDaughters) {
            thrust::copy_if(
                fWeights.begin(), fWeights.end(), fAccRejFlags.begin(), _Events->fWeights.begin(), isAccepted());

            thrust::copy_if(fAccRejFlags.begin(),
                            fAccRejFlags.end(),
                            fAccRejFlags.begin(),
                            _Events->fAccRejFlags.begin(),
                            isAccepted());
        }
    }

#else

    cudaStream_t s[fNDaughters + 1];

    for(GInt_t d = 0; d <= fNDaughters; d++) {
        CUDA_CHECK_RETURN(cudaStreamCreateWithFlags(&s[d], cudaStreamNonBlocking));
    }

    thrust::copy_if(thrust::cuda::par.on(s[fNDaughters]),
                    fWeights.begin(),
                    fWeights.end(),
                    fAccRejFlags.begin(),
                    _Events->fWeights.begin(),
                    isAccepted());

    thrust::copy_if(thrust::cuda::par.on(s[fNDaughters]),
                    fAccRejFlags.begin(),
                    fAccRejFlags.end(),
                    fAccRejFlags.begin(),
                    _Events->fAccRejFlags.begin(),
                    isAccepted());

    for(GInt_t d = 0; d < fNDaughters; d++) {
        thrust::copy_if(thrust::cuda::par.on(s[d]),
                        fDaughters[d].begin(),
                        fDaughters[d].end(),
                        fAccRejFlags.begin(),
                        _Events->fDaughters[d].begin(),
                        isAccepted());
    }

    cudaDeviceSynchronize();
    for(GInt_t d = 0; d <= fNDaughters; d++)
        cudaStreamDestroy(s[d]);

#endif
}

inline void PhaseSpace::Export(Events *_Events) {
    /**
     * Export the events and all related information to an Events object properly initialized.
     */
    _Events->fMaxWeight = fMaxWeight;

#if MCBOOSTER_BACKEND != CUDA
#if MCBOOSTER_BACKEND != CPP
#pragma omp parallel num_threads(fNDaughters + 1)
    {
#else
    for(size_t val = 0; val < fNDaughters; val++) {
        auto omp_get_thread_num = [&val]() { return val; };
#endif

        if(omp_get_thread_num() < fNDaughters) {
            thrust::copy(fDaughters[omp_get_thread_num()].begin(),
                         fDaughters[omp_get_thread_num()].end(),
                         _Events->fDaughters[omp_get_thread_num()].begin());
        }

        if(omp_get_thread_num() == fNDaughters) {
            thrust::copy(fWeights.begin(), fWeights.end(), _Events->fWeights.begin());

            thrust::copy(fAccRejFlags.begin(), fAccRejFlags.end(), _Events->fAccRejFlags.begin());
        }
    }

#else

    cudaStream_t s[fNDaughters + 1];

    for(GInt_t d = 0; d <= fNDaughters; d++) {
        CUDA_CHECK_RETURN(cudaStreamCreateWithFlags(&s[d], cudaStreamNonBlocking));
    }

    cudaMemcpyAsync(thrust::raw_pointer_cast(_Events->fWeights.data()),
                    thrust::raw_pointer_cast(fWeights.data()),
                    fWeights.size() * sizeof(GReal_t),
                    cudaMemcpyDeviceToHost,
                    s[fNDaughters]);

    cudaMemcpyAsync(thrust::raw_pointer_cast(_Events->fAccRejFlags.data()),
                    thrust::raw_pointer_cast(fAccRejFlags.data()),
                    fAccRejFlags.size() * sizeof(GBool_t),
                    cudaMemcpyDeviceToHost,
                    s[fNDaughters]);

    for(GInt_t d = 0; d < fNDaughters; d++) {
        cudaMemcpyAsync(thrust::raw_pointer_cast(_Events->fDaughters[d].data()),
                        thrust::raw_pointer_cast(fDaughters[d].data()),
                        fDaughters[d].size() * sizeof(Vector4R),
                        cudaMemcpyDeviceToHost,
                        s[d]);
    }

    cudaDeviceSynchronize();
    for(GInt_t d = 0; d <= fNDaughters; d++)
        cudaStreamDestroy(s[d]);

#endif
}

inline void PhaseSpace::Generate(const Vector4R fMother) {
/**
 * Run the generator and calculate the maximum weight. It takes as input the fourvector of the mother particle
 * in any system of reference. The daughters will be generated in this system.
 */

#if MCBOOSTER_BACKEND == CUDA
    cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
#endif
    /* random number generation */

    RND_Time = 0;
    // create iterators
    thrust::counting_iterator<GLong_t> first(0);
    thrust::counting_iterator<GLong_t> last = first + fNEvents;

    // Vai!!!

    /* event generation */
    timespec time_event_start, time_event_end;
    clock_gettime(TIMER, &time_event_start);

    switch(fNDaughters) {
    case 2:

        thrust::transform(first,
                          last,
                          thrust::make_zip_iterator(thrust::make_tuple(fDaughters[0].begin(), fDaughters[1].begin())),
                          fWeights.begin(),
                          DecayMother(fMother, fMasses, fNDaughters, fSeed));

        break;

    case 3:

        thrust::transform(first,
                          last,
                          thrust::make_zip_iterator(
                              thrust::make_tuple(fDaughters[0].begin(), fDaughters[1].begin(), fDaughters[2].begin())),
                          fWeights.begin(),
                          DecayMother(fMother, fMasses, fNDaughters, fSeed));

        break;
    case 4:

        thrust::transform(
            first,
            last,
            thrust::make_zip_iterator(thrust::make_tuple(
                fDaughters[0].begin(), fDaughters[1].begin(), fDaughters[2].begin(), fDaughters[3].begin())),
            fWeights.begin(),
            DecayMother(fMother, fMasses, fNDaughters, fSeed));

        break;
    case 5:

        thrust::transform(first,
                          last,
                          thrust::make_zip_iterator(thrust::make_tuple(fDaughters[0].begin(),
                                                                       fDaughters[1].begin(),
                                                                       fDaughters[2].begin(),
                                                                       fDaughters[3].begin(),
                                                                       fDaughters[4].begin())),
                          fWeights.begin(),
                          DecayMother(fMother, fMasses, fNDaughters, fSeed));

        //}
        break;
    case 6:

        thrust::transform(first,
                          last,
                          thrust::make_zip_iterator(thrust::make_tuple(fDaughters[0].begin(),
                                                                       fDaughters[1].begin(),
                                                                       fDaughters[2].begin(),
                                                                       fDaughters[3].begin(),
                                                                       fDaughters[4].begin(),
                                                                       fDaughters[5].begin())),
                          fWeights.begin(),
                          DecayMother(fMother, fMasses, fNDaughters, fSeed));

        break;
    case 7:

        thrust::transform(first,
                          last,
                          thrust::make_zip_iterator(thrust::make_tuple(fDaughters[0].begin(),
                                                                       fDaughters[1].begin(),
                                                                       fDaughters[2].begin(),
                                                                       fDaughters[3].begin(),
                                                                       fDaughters[4].begin(),
                                                                       fDaughters[5].begin(),
                                                                       fDaughters[6].begin())),
                          fWeights.begin(),
                          DecayMother(fMother, fMasses, fNDaughters, fSeed));

        break;
    case 8:

        thrust::transform(first,
                          last,
                          thrust::make_zip_iterator(thrust::make_tuple(fDaughters[0].begin(),
                                                                       fDaughters[1].begin(),
                                                                       fDaughters[2].begin(),
                                                                       fDaughters[3].begin(),
                                                                       fDaughters[4].begin(),
                                                                       fDaughters[5].begin(),
                                                                       fDaughters[6].begin(),
                                                                       fDaughters[7].begin())),
                          fWeights.begin(),
                          DecayMother(fMother, fMasses, fNDaughters, fSeed));

        break;
    case 9:

        thrust::transform(first,
                          last,
                          thrust::make_zip_iterator(thrust::make_tuple(fDaughters[0].begin(),
                                                                       fDaughters[1].begin(),
                                                                       fDaughters[2].begin(),
                                                                       fDaughters[3].begin(),
                                                                       fDaughters[4].begin(),
                                                                       fDaughters[5].begin(),
                                                                       fDaughters[6].begin(),
                                                                       fDaughters[7].begin(),
                                                                       fDaughters[8].begin())),
                          fWeights.begin(),
                          DecayMother(fMother, fMasses, fNDaughters, fSeed));

        break;
    }

    clock_gettime(TIMER, &time_event_end);
    EVT_Time = ((GReal_t)(time_diff(time_event_start, time_event_end).tv_sec
                          + time_diff(time_event_start, time_event_end).tv_nsec * 1.0e-9));

    // setting maximum weight
    RealVector_d::iterator w = thrust::max_element(fWeights.begin(), fWeights.end());
    fMaxWeight               = *w;
}

inline void PhaseSpace::Generate(Particles_d fMothers) {
/**
 * Run the generator and calculate the maximum weight. It takes as input the device vector with the four-vectors of the
 * mother particle
 * in any system of reference. The daughters will be generated in this system.
 */

#if MCBOOSTER_BACKEND == CUDA
    cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
#endif

    if(fNEvents < fMothers.size())
        std::cout << "fNEvents != fMothers.size()" << std::endl;

    /* random number generation */
    /*
    timespec time_rnd_start, time_rnd_end;
    clock_gettime(TIMER, &time_rnd_start);


    clock_gettime(TIMER, &time_rnd_end);

    RND_Time = ((GReal_t) (time_diff(time_rnd_start, time_rnd_end).tv_sec
            + time_diff(time_rnd_start, time_rnd_end).tv_nsec * 1.0e-9));
    */

    /* event generation */
    timespec time_event_start, time_event_end;
    clock_gettime(TIMER, &time_event_start);

    RND_Time = 0.0;
    // create iterators
    thrust::counting_iterator<GLong_t> first(0);
    thrust::counting_iterator<GLong_t> last = first + fNEvents;

    switch(fNDaughters) {
    case 2:

        thrust::transform(first,
                          last,
                          thrust::make_zip_iterator(
                              thrust::make_tuple(fMothers.begin(), fDaughters[0].begin(), fDaughters[1].begin())),
                          fWeights.begin(),
                          DecayMothers(fMasses, fNDaughters, fSeed));

        break;

    case 3:

        thrust::transform(first,
                          last,
                          thrust::make_zip_iterator(thrust::make_tuple(
                              fMothers.begin(), fDaughters[0].begin(), fDaughters[1].begin(), fDaughters[2].begin())),
                          fWeights.begin(),
                          DecayMothers(fMasses, fNDaughters, fSeed));

        break;
    case 4:

        thrust::transform(first,
                          last,
                          thrust::make_zip_iterator(thrust::make_tuple(fMothers.begin(),
                                                                       fDaughters[0].begin(),
                                                                       fDaughters[1].begin(),
                                                                       fDaughters[2].begin(),
                                                                       fDaughters[3].begin())),
                          fWeights.begin(),
                          DecayMothers(fMasses, fNDaughters, fSeed));

        break;
    case 5:

        thrust::transform(first,
                          last,
                          thrust::make_zip_iterator(thrust::make_tuple(fMothers.begin(),
                                                                       fDaughters[0].begin(),
                                                                       fDaughters[1].begin(),
                                                                       fDaughters[2].begin(),
                                                                       fDaughters[3].begin(),
                                                                       fDaughters[4].begin())),
                          fWeights.begin(),
                          DecayMothers(fMasses, fNDaughters, fSeed));

        break;
    case 6:

        thrust::transform(first,
                          last,
                          thrust::make_zip_iterator(thrust::make_tuple(fMothers.begin(),
                                                                       fDaughters[0].begin(),
                                                                       fDaughters[1].begin(),
                                                                       fDaughters[2].begin(),
                                                                       fDaughters[3].begin(),
                                                                       fDaughters[4].begin(),
                                                                       fDaughters[5].begin())),
                          fWeights.begin(),
                          DecayMothers(fMasses, fNDaughters, fSeed));

        break;
    case 7:

        thrust::transform(first,
                          last,
                          thrust::make_zip_iterator(thrust::make_tuple(fMothers.begin(),
                                                                       fDaughters[0].begin(),
                                                                       fDaughters[1].begin(),
                                                                       fDaughters[2].begin(),
                                                                       fDaughters[3].begin(),
                                                                       fDaughters[4].begin(),
                                                                       fDaughters[5].begin(),
                                                                       fDaughters[6].begin())),
                          fWeights.begin(),
                          DecayMothers(fMasses, fNDaughters, fSeed));

        break;
    case 8:

        thrust::transform(first,
                          last,
                          thrust::make_zip_iterator(thrust::make_tuple(fMothers.begin(),
                                                                       fDaughters[0].begin(),
                                                                       fDaughters[1].begin(),
                                                                       fDaughters[2].begin(),
                                                                       fDaughters[3].begin(),
                                                                       fDaughters[4].begin(),
                                                                       fDaughters[5].begin(),
                                                                       fDaughters[6].begin(),
                                                                       fDaughters[7].begin())),
                          fWeights.begin(),
                          DecayMothers(fMasses, fNDaughters, fSeed));

        break;
    case 9:

        thrust::transform(first,
                          last,
                          thrust::make_zip_iterator(thrust::make_tuple(fMothers.begin(),
                                                                       fDaughters[0].begin(),
                                                                       fDaughters[1].begin(),
                                                                       fDaughters[2].begin(),
                                                                       fDaughters[3].begin(),
                                                                       fDaughters[4].begin(),
                                                                       fDaughters[5].begin(),
                                                                       fDaughters[6].begin(),
                                                                       fDaughters[7].begin(),
                                                                       fDaughters[8].begin())),
                          fWeights.begin(),
                          DecayMothers(fMasses, fNDaughters, fSeed));

        break;
    }

    clock_gettime(TIMER, &time_event_end);
    EVT_Time = ((GReal_t)(time_diff(time_event_start, time_event_end).tv_sec
                          + time_diff(time_event_start, time_event_end).tv_nsec * 1.0e-9));

    // setting maximum weight
    RealVector_d::iterator w = thrust::max_element(fWeights.begin(), fWeights.end());
    fMaxWeight               = *w;
}
}
#endif /* GENERATE_H_ */
