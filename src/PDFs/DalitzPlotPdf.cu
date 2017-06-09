#include "goofit/PDFs/physics/DalitzPlotPdf.h"
#include "goofit/Error.h"

#include <thrust/complex.h>
#include <thrust/transform_reduce.h>

namespace GooFit {

const int resonanceOffset_DP = 4; // Offset of the first resonance into the parameter index array
// Offset is number of parameters, constant index, number of resonances (not calculable
// from nP because we don't know what the efficiency might need), and cache index. Efficiency
// parameters are after the resonance information.

// The function of this array is to hold all the cached waves; specific
// waves are recalculated when the corresponding resonance mass or width
// changes. Note that in a multithread environment each thread needs its
// own cache, hence the '10'. Ten threads should be enough for anyone!

// NOTE: This is does not support ten instances (ten threads) of resoncances now, only one set of resonances.
__device__ thrust::complex<fptype> *cResonances[16];

__device__ inline int parIndexFromResIndex_DP(int resIndex) { return resonanceOffset_DP + resIndex * resonanceSize; }

__device__ thrust::complex<fptype>
device_DalitzPlot_calcIntegrals(fptype m12, fptype m13, int res_i, int res_j, fptype *p, unsigned int *indices) {
    // Calculates BW_i(m12, m13) * BW_j^*(m12, m13).
    // This calculation is in a separate function so
    // it can be cached. Note that this function expects
    // to be called on a normalisation grid, not on
    // observed points, that's why it doesn't use
    // cResonances. No need to cache the values at individual
    // grid points - we only care about totals.
    fptype motherMass = RO_CACHE(functorConstants[RO_CACHE(indices[1]) + 0]);
    fptype daug1Mass  = RO_CACHE(functorConstants[RO_CACHE(indices[1]) + 1]);
    fptype daug2Mass  = RO_CACHE(functorConstants[RO_CACHE(indices[1]) + 2]);
    fptype daug3Mass  = RO_CACHE(functorConstants[RO_CACHE(indices[1]) + 3]);

    thrust::complex<fptype> ret;

    if(!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass))
        return ret;

    fptype m23
        = motherMass * motherMass + daug1Mass * daug1Mass + daug2Mass * daug2Mass + daug3Mass * daug3Mass - m12 - m13;

    int parameter_i       = parIndexFromResIndex_DP(res_i);
    unsigned int functn_i = RO_CACHE(indices[parameter_i + 2]);
    unsigned int params_i = RO_CACHE(indices[parameter_i + 3]);
    ret                   = getResonanceAmplitude(m12, m13, m23, functn_i, params_i);

    int parameter_j       = parIndexFromResIndex_DP(res_j);
    unsigned int functn_j = RO_CACHE(indices[parameter_j + 2]);
    unsigned int params_j = RO_CACHE(indices[parameter_j + 3]);
    ret *= conj(getResonanceAmplitude(m12, m13, m23, functn_j, params_j));

    return ret;
}

__device__ fptype device_DalitzPlot(fptype *evt, fptype *p, unsigned int *indices) {
    fptype motherMass = RO_CACHE(functorConstants[RO_CACHE(indices[1]) + 0]);
    fptype daug1Mass  = RO_CACHE(functorConstants[RO_CACHE(indices[1]) + 1]);
    fptype daug2Mass  = RO_CACHE(functorConstants[RO_CACHE(indices[1]) + 2]);
    fptype daug3Mass  = RO_CACHE(functorConstants[RO_CACHE(indices[1]) + 3]);

    fptype m12 = RO_CACHE(evt[RO_CACHE(indices[2 + RO_CACHE(indices[0])])]);
    fptype m13 = RO_CACHE(evt[RO_CACHE(indices[3 + RO_CACHE(indices[0])])]);

    if(!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass))
        return 0;

    fptype evtIndex = RO_CACHE(evt[RO_CACHE(indices[4 + RO_CACHE(indices[0])])]);

    auto evtNum = static_cast<int>(floor(0.5 + evtIndex));

    thrust::complex<fptype> totalAmp(0, 0);
    unsigned int numResonances = RO_CACHE(indices[2]);
    // unsigned int cacheToUse    = RO_CACHE(indices[3]);

    for(int i = 0; i < numResonances; ++i) {
        int paramIndex              = parIndexFromResIndex_DP(i);
        thrust::complex<fptype> amp = thrust::complex<fptype>(RO_CACHE(p[RO_CACHE(indices[paramIndex + 0])]),
                                                              RO_CACHE(p[RO_CACHE(indices[paramIndex + 1])]));

        // potential performance improvement by
        // double2 me = RO_CACHE(reinterpret_cast<double2*> (cResonances[i][evtNum]));
        thrust::complex<fptype> me = RO_CACHE(cResonances[i][evtNum]);

        // thrust::complex<fptype> matrixelement((cResonances[cacheToUse][evtNum*numResonances + i]).real,
        //				     (cResonances[cacheToUse][evtNum*numResonances + i]).imag);
        // thrust::complex<fptype> matrixelement (me[0], me[1]);

        totalAmp += amp * me;
    }

    fptype ret         = thrust::norm(totalAmp);
    int effFunctionIdx = parIndexFromResIndex_DP(numResonances);
    fptype eff         = callFunction(evt, RO_CACHE(indices[effFunctionIdx]), RO_CACHE(indices[effFunctionIdx + 1]));
    ret *= eff;

    // printf("DalitzPlot evt %i zero: %i %i %f (%f, %f).\n", evtNum, numResonances, effFunctionIdx, eff, totalAmp.real,
    // totalAmp.imag);

    return ret;
}

__device__ device_function_ptr ptr_to_DalitzPlot = device_DalitzPlot;

__host__ DalitzPlotPdf::DalitzPlotPdf(
    std::string n, Variable *m12, Variable *m13, CountingVariable *eventNumber, DecayInfo *decay, GooPdf *efficiency)
    : GooPdf(nullptr, n)
    , decayInfo(decay)
    , _m12(m12)
    , _m13(m13)
    , dalitzNormRange(nullptr)
    //, cachedWaves(0)
    , integrals(nullptr)
    , forceRedoIntegrals(true)
    , totalEventSize(3) // Default 3 = m12, m13, evtNum
    , cacheToUse(0)
    , integrators(nullptr)
    , calculators(nullptr) {
    registerObservable(_m12);
    registerObservable(_m13);
    registerObservable(eventNumber);

    fptype decayConstants[5];

    for(auto &cachedWave : cachedWaves)
        cachedWave = nullptr;

    std::vector<unsigned int> pindices;
    pindices.push_back(registerConstants(5));
    decayConstants[0] = decayInfo->motherMass;
    decayConstants[1] = decayInfo->daug1Mass;
    decayConstants[2] = decayInfo->daug2Mass;
    decayConstants[3] = decayInfo->daug3Mass;
    decayConstants[4] = decayInfo->meson_radius;
    MEMCPY_TO_SYMBOL(
        functorConstants, decayConstants, 5 * sizeof(fptype), cIndex * sizeof(fptype), cudaMemcpyHostToDevice);

    pindices.push_back(decayInfo->resonances.size());
    static int cacheCount = 0;
    cacheToUse            = cacheCount++;
    pindices.push_back(cacheToUse);

    for(auto &resonance : decayInfo->resonances) {
        pindices.push_back(registerParameter(resonance->amp_real));
        pindices.push_back(registerParameter(resonance->amp_imag));
        pindices.push_back(resonance->getFunctionIndex());
        pindices.push_back(resonance->getParameterIndex());
        resonance->setConstantIndex(cIndex);
        components.push_back(resonance);
    }

    pindices.push_back(efficiency->getFunctionIndex());
    pindices.push_back(efficiency->getParameterIndex());
    components.push_back(efficiency);

    GET_FUNCTION_ADDR(ptr_to_DalitzPlot);
    initialize(pindices);

    redoIntegral = new bool[decayInfo->resonances.size()];
    cachedMasses = new fptype[decayInfo->resonances.size()];
    cachedWidths = new fptype[decayInfo->resonances.size()];
    integrals    = new thrust::complex<fptype> **[decayInfo->resonances.size()];
    integrators  = new SpecialResonanceIntegrator **[decayInfo->resonances.size()];
    calculators  = new SpecialResonanceCalculator *[decayInfo->resonances.size()];

    for(int i = 0; i < decayInfo->resonances.size(); ++i) {
        redoIntegral[i] = true;
        cachedMasses[i] = -1;
        cachedWidths[i] = -1;
        integrators[i]  = new SpecialResonanceIntegrator *[decayInfo->resonances.size()];
        calculators[i]  = new SpecialResonanceCalculator(parameters, i);
        integrals[i]    = new thrust::complex<fptype> *[decayInfo->resonances.size()];

        for(int j = 0; j < decayInfo->resonances.size(); ++j) {
            integrals[i][j]   = new thrust::complex<fptype>(0, 0);
            integrators[i][j] = new SpecialResonanceIntegrator(parameters, i, j);
        }
    }

    addSpecialMask(PdfBase::ForceSeparateNorm);
}

__host__ void DalitzPlotPdf::setDataSize(unsigned int dataSize, unsigned int evtSize) {
    // Default 3 is m12, m13, evtNum
    totalEventSize = evtSize;
    if(totalEventSize < 3)
        throw GooFit::GeneralError("totalEventSize {} must be 3 or more", totalEventSize);

    // if (cachedWaves) delete cachedWaves;
    if(cachedWaves[0]) {
        for(auto &cachedWave : cachedWaves)
            delete cachedWave;
    }

    numEntries = dataSize;

    for(int i = 0; i < 16; i++) {
#ifdef GOOFIT_MPI
        cachedWaves[i] = new thrust::device_vector<thrust::complex<fptype>>(m_iEventsPerTask);
#else
        cachedWaves[i] = new thrust::device_vector<thrust::complex<fptype>>(dataSize);
#endif
        void *dummy = thrust::raw_pointer_cast(cachedWaves[i]->data());
        MEMCPY_TO_SYMBOL(cResonances,
                         &dummy,
                         sizeof(thrust::complex<fptype> *),
                         i * sizeof(thrust::complex<fptype> *),
                         cudaMemcpyHostToDevice);
    }

    setForceIntegrals();
}

__host__ fptype DalitzPlotPdf::normalize() const {
    recursiveSetNormalisation(1); // Not going to normalize efficiency,
    // so set normalisation factor to 1 so it doesn't get multiplied by zero.
    // Copy at this time to ensure that the SpecialResonanceCalculators, which need the efficiency,
    // don't get zeroes through multiplying by the normFactor.
    MEMCPY_TO_SYMBOL(normalisationFactors, host_normalisation, totalParams * sizeof(fptype), 0, cudaMemcpyHostToDevice);

    int totalBins = _m12->getNumBins() * _m13->getNumBins();

    if(!dalitzNormRange) {
        gooMalloc((void **)&dalitzNormRange, 6 * sizeof(fptype));

        auto *host_norms = new fptype[6];
        host_norms[0]    = _m12->getLowerLimit();
        host_norms[1]    = _m12->getUpperLimit();
        host_norms[2]    = _m12->getNumBins();
        host_norms[3]    = _m13->getLowerLimit();
        host_norms[4]    = _m13->getUpperLimit();
        host_norms[5]    = _m13->getNumBins();
        MEMCPY(dalitzNormRange, host_norms, 6 * sizeof(fptype), cudaMemcpyHostToDevice);
        delete[] host_norms;
    }

    for(unsigned int i = 0; i < decayInfo->resonances.size(); ++i) {
        redoIntegral[i] = forceRedoIntegrals;

        if(!(decayInfo->resonances[i]->parametersChanged()))
            continue;

        redoIntegral[i] = true;
    }

    forceRedoIntegrals = false;

    // Only do this bit if masses or widths have changed.
    thrust::constant_iterator<fptype *> arrayAddress(dalitzNormRange);
    thrust::counting_iterator<int> binIndex(0);

    // NB, SpecialResonanceCalculator assumes that fit is unbinned!
    // And it needs to know the total event size, not just observables
    // for this particular PDF component.
    thrust::constant_iterator<fptype *> dataArray(dev_event_array);
    thrust::constant_iterator<int> eventSize(totalEventSize);
    thrust::counting_iterator<int> eventIndex(0);

    for(int i = 0; i < decayInfo->resonances.size(); ++i) {
        if(redoIntegral[i]) {
#ifdef GOOFIT_MPI
            thrust::transform(
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex + m_iEventsPerTask, arrayAddress, eventSize)),
                strided_range<thrust::device_vector<thrust::complex<fptype>>::iterator>(
                    cachedWaves[i]->begin(), cachedWaves[i]->end(), 1)
                    .begin(),
                *(calculators[i]));
#else
            thrust::transform(
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, arrayAddress, eventSize)),
                strided_range<thrust::device_vector<thrust::complex<fptype>>::iterator>(
                    cachedWaves[i]->begin(), cachedWaves[i]->end(), 1)
                    .begin(),
                *(calculators[i]));
#endif
        }

        // Possibly this can be done more efficiently by exploiting symmetry?
        for(int j = 0; j < decayInfo->resonances.size(); ++j) {
            if((!redoIntegral[i]) && (!redoIntegral[j]))
                continue;

            thrust::complex<fptype> dummy(0, 0);
            thrust::plus<thrust::complex<fptype>> complexSum;
            (*(integrals[i][j])) = thrust::transform_reduce(
                thrust::make_zip_iterator(thrust::make_tuple(binIndex, arrayAddress)),
                thrust::make_zip_iterator(thrust::make_tuple(binIndex + totalBins, arrayAddress)),
                *(integrators[i][j]),
                dummy,
                complexSum);
        }
    }

    // End of time-consuming integrals.
    thrust::complex<fptype> sumIntegral(0, 0);

    for(unsigned int i = 0; i < decayInfo->resonances.size(); ++i) {
        int param_i = parameters + resonanceOffset_DP + resonanceSize * i;
        thrust::complex<fptype> amplitude_i(host_params[host_indices[param_i]], host_params[host_indices[param_i + 1]]);

        for(unsigned int j = 0; j < decayInfo->resonances.size(); ++j) {
            int param_j = parameters + resonanceOffset_DP + resonanceSize * j;
            thrust::complex<fptype> amplitude_j(host_params[host_indices[param_j]],
                                                -host_params[host_indices[param_j + 1]]);
            // Notice complex conjugation
            // printf("%f %f %f %f %f %f\n", amplitude_i.real(), amplitude_i.imag(), amplitude_j.real(),
            // amplitude_j.imag(), (*(integrals[i][j])).real, (*(integrals[i][j])).imag );
            sumIntegral += amplitude_i * amplitude_j * (*(integrals[i][j]));
        }
    }

    fptype ret           = sumIntegral.real(); // That complex number is a square, so it's fully real
    double binSizeFactor = 1;
    binSizeFactor *= _m12->getBinSize();
    binSizeFactor *= _m13->getBinSize();
    ret *= binSizeFactor;

    host_normalisation[parameters] = 1.0 / ret;
    // printf("%f %f\n", ret, binSizeFactor);
    return ret;
}

SpecialResonanceIntegrator::SpecialResonanceIntegrator(int pIdx, unsigned int ri, unsigned int rj)
    : resonance_i(ri)
    , resonance_j(rj)
    , parameters(pIdx) {}

__device__ thrust::complex<fptype> SpecialResonanceIntegrator::operator()(thrust::tuple<int, fptype *> t) const {
    // Bin index, base address [lower, upper,getNumBins]
    // Notice that this is basically MetricTaker::operator (binned) with the special-case knowledge
    // that event size is two, and that the function to call is dev_DalitzPlot_calcIntegrals.

    int globalBinNumber  = thrust::get<0>(t);
    fptype lowerBoundM12 = thrust::get<1>(t)[0];
    fptype upperBoundM12 = thrust::get<1>(t)[1];
    auto numBinsM12      = static_cast<int>(floor(thrust::get<1>(t)[2] + 0.5));
    int binNumberM12     = globalBinNumber % numBinsM12;
    fptype binCenterM12  = upperBoundM12 - lowerBoundM12;
    binCenterM12 /= numBinsM12;
    binCenterM12 *= (binNumberM12 + 0.5);
    binCenterM12 += lowerBoundM12;

    globalBinNumber /= numBinsM12;
    fptype lowerBoundM13 = thrust::get<1>(t)[3];
    fptype upperBoundM13 = thrust::get<1>(t)[4];
    auto numBinsM13      = static_cast<int>(floor(thrust::get<1>(t)[5] + 0.5));
    fptype binCenterM13  = upperBoundM13 - lowerBoundM13;
    binCenterM13 /= numBinsM13;
    binCenterM13 *= (globalBinNumber + 0.5);
    binCenterM13 += lowerBoundM13;

    unsigned int *indices = paramIndices + parameters;
    thrust::complex<fptype> ret
        = device_DalitzPlot_calcIntegrals(binCenterM12, binCenterM13, resonance_i, resonance_j, cudaArray, indices);

    fptype fakeEvt[10]; // Need room for many observables in case m12 or m13 were assigned a high index in an
                        // event-weighted fit.
    fakeEvt[indices[indices[0] + 2 + 0]] = binCenterM12;
    fakeEvt[indices[indices[0] + 2 + 1]] = binCenterM13;
    unsigned int numResonances           = indices[2];
    int effFunctionIdx                   = parIndexFromResIndex_DP(numResonances);
    fptype eff                           = callFunction(fakeEvt, indices[effFunctionIdx], indices[effFunctionIdx + 1]);

    // Multiplication by eff, not sqrt(eff), is correct:
    // These complex numbers will not be squared when they
    // go into the integrals. They've been squared already,
    // as it were.
    ret *= eff;
    // printf("ret %f %f %f %f %f\n",binCenterM12, binCenterM13, ret.real, ret.imag, eff );
    return ret;
}

SpecialResonanceCalculator::SpecialResonanceCalculator(int pIdx, unsigned int res_idx)
    : resonance_i(res_idx)
    , parameters(pIdx) {}

__device__ thrust::complex<fptype> SpecialResonanceCalculator::operator()(thrust::tuple<int, fptype *, int> t) const {
    // Calculates the BW values for a specific resonance.
    thrust::complex<fptype> ret;
    int evtNum  = thrust::get<0>(t);
    fptype *evt = thrust::get<1>(t) + (evtNum * thrust::get<2>(t));

    unsigned int *indices = paramIndices + parameters; // Jump to DALITZPLOT position within parameters array
    fptype m12            = evt[indices[2 + indices[0]]];
    fptype m13            = evt[indices[3 + indices[0]]];

    fptype motherMass = functorConstants[indices[1] + 0];
    fptype daug1Mass  = functorConstants[indices[1] + 1];
    fptype daug2Mass  = functorConstants[indices[1] + 2];
    fptype daug3Mass  = functorConstants[indices[1] + 3];

    if(!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass))
        return ret;

    fptype m23
        = motherMass * motherMass + daug1Mass * daug1Mass + daug2Mass * daug2Mass + daug3Mass * daug3Mass - m12 - m13;

    int parameter_i
        = parIndexFromResIndex_DP(resonance_i); // Find position of this resonance relative to DALITZPLOT start

    unsigned int functn_i = indices[parameter_i + 2];
    unsigned int params_i = indices[parameter_i + 3];

    ret = getResonanceAmplitude(m12, m13, m23, functn_i, params_i);
    // printf("Amplitude %f %f %f (%f, %f)\n ", m12, m13, m23, ret.real, ret.imag);
    return ret;
}

} // namespace GooFit
