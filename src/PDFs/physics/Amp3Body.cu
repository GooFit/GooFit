#include <goofit/Error.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/Amp3Body.h>
#include <goofit/PDFs/physics/detail/SpecialResonanceCalculator.h>
#include <goofit/PDFs/physics/detail/SpecialResonanceIntegrator.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>
#include <goofit/detail/Complex.h>

#include <thrust/copy.h>
#include <thrust/transform_reduce.h>

#include <array>

namespace GooFit {

// Functor used for fit fraction sum
struct CoefSumFunctor {
    fpcomplex coef_i;
    fpcomplex coef_j;

    CoefSumFunctor(fpcomplex coef_i, fpcomplex coef_j)
        : coef_i(coef_i)
        , coef_j(coef_j) {}

    __device__ fptype operator()(thrust::tuple<fpcomplex, fpcomplex> val) {
        return (coef_i * thrust::conj<fptype>(coef_j) * thrust::get<0>(val) * thrust::conj<fptype>(thrust::get<1>(val)))
            .real();
    }
};

constexpr int resonanceOffset_DP = 4; // Offset of the first resonance into the parameter index array
// Offset is number of parameters, constant index, number of resonances (not calculable
// from nP because we don't know what the efficiency might need), and cache index. Efficiency
// parameters are after the resonance information.

// The function of this array is to hold all the cached waves; specific
// waves are recalculated when the corresponding resonance mass or width
// changes. Note that in a multithread environment each thread needs its
// own cache, hence the '10'. Ten threads should be enough for anyone!

// NOTE: This is does not support ten instances (ten threads) of resoncances now, only one set of resonances.
__device__ fpcomplex *cResonances[16];

__device__ inline int parIndexFromResIndex_DP(int resIndex) { return resonanceOffset_DP + resIndex * resonanceSize; }

__device__ fptype device_DalitzPlot(fptype *evt, ParameterContainer &pc) {
    int num_obs = pc.getNumObservables();
    int id_m12  = pc.getObservable(0);
    int id_m13  = pc.getObservable(1);
    int id_num  = pc.getObservable(2);

    fptype m12 = RO_CACHE(evt[id_m12]);
    fptype m13 = RO_CACHE(evt[id_m13]);

    unsigned int numResonances = pc.getConstant(0);
    // unsigned int cacheToUse    = pc.getConstant(1);

    if(!inDalitz(m12, m13, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass)) {
        pc.incrementIndex(1, numResonances * 2, 2, num_obs, 1);

        // loop over resonances and efficiency functions
        for(int i = 0; i < numResonances; i++)
            pc.incrementIndex();

        // increment the efficiency function
        pc.incrementIndex();
        return 0;
    }

    fptype evtIndex = RO_CACHE(evt[id_num]);

    auto evtNum = static_cast<int>(floor(0.5 + evtIndex));

    fpcomplex totalAmp(0, 0);

    for(int i = 0; i < numResonances; ++i) {
        fpcomplex amp = fpcomplex(pc.getParameter(i * 2), pc.getParameter(i * 2 + 1));

        // potential performance improvement by
        // double2 *t = RO_CACHE(reinterpret_cast<double2*> (&(cResonances[i][evtNum])));
        // fpcomplex me(t->x, t->y);
        // fpcomplex me = RO_CACHE(cResonances[i][evtNum]);

        // double2 *ptr = reinterpret_cast<double2*> (cResonances[i][evtNum]);

        // double2 v = RO_CACHE(resPtr[evtNum]);

        // fptype me_real = cResonances[i][evtNum].real();
        // fptype me_imag = cResonances[i][evtNum].imag();
        // fpcomplex me = cResonances[i][evtNum];
        // fpcomplex me (me_real, me_imag);
        fpcomplex me = RO_CACHE(cResonances[i][evtNum]);

        totalAmp += amp * me;
    }

    fptype ret = thrust::norm(totalAmp);
    pc.incrementIndex(1, numResonances * 2, 2, num_obs, 1);

    // loop to efficiency idx
    for(int i = 0; i < numResonances; i++)
        pc.incrementIndex();

    fptype eff = callFunction(evt, pc);
    ret *= eff;

    return ret;
}

__device__ device_function_ptr ptr_to_DalitzPlot = device_DalitzPlot;

__host__ Amp3Body::Amp3Body(
    std::string n, Observable m12, Observable m13, EventNumber eventNumber, DecayInfo3 decay, GooPdf *efficiency)
    : Amp3BodyBase("Amp3Body", n, m12, m13, eventNumber)
    , decayInfo(decay)
    , _m12(m12)
    , _m13(m13)
    , _eventNumber(eventNumber)
    , dalitzNormRange(nullptr)
    //, cachedWaves(0)
    , integrals(nullptr)
    , forceRedoIntegrals(true)
    , totalEventSize(3) // Default 3 = m12, m13, evtNum
    , cacheToUse(0)
    , integrators(nullptr)
    , calculators(nullptr) {
    for(auto &cachedWave : cachedWaves)
        cachedWave = nullptr;

    // Passing values to the defined constants.  Rather than push into list, which means each resonance
    MEMCPY_TO_SYMBOL(c_motherMass, &decay.motherMass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug1Mass, &decay.daug1Mass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug2Mass, &decay.daug2Mass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug3Mass, &decay.daug3Mass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_meson_radius, &decay.meson_radius, sizeof(fptype), 0, cudaMemcpyHostToDevice);

    // registered to 0 position
    registerConstant(decayInfo.resonances.size());
    static int cacheCount = 0;
    cacheToUse            = cacheCount++;
    // registered to 1 position
    registerConstant(cacheToUse);

    for(auto &resonance : decayInfo.resonances) {
        // registering 2 parameters
        registerParameter(resonance->amp_real);
        registerParameter(resonance->amp_imag);
        components.push_back(resonance);
    }

    components.push_back(efficiency);

    registerFunction("ptr_to_DalitzPlot", ptr_to_DalitzPlot);

    initialize();

    redoIntegral = new bool[decayInfo.resonances.size()];
    cachedMasses = new fptype[decayInfo.resonances.size()];
    cachedWidths = new fptype[decayInfo.resonances.size()];
    integrals    = new fpcomplex **[decayInfo.resonances.size()];
    integrators  = new SpecialResonanceIntegrator **[decayInfo.resonances.size()];
    calculators  = new SpecialResonanceCalculator *[decayInfo.resonances.size()];

    for(int i = 0; i < decayInfo.resonances.size(); ++i) {
        redoIntegral[i] = true;
        cachedMasses[i] = -1;
        cachedWidths[i] = -1;
        integrators[i]  = new SpecialResonanceIntegrator *[decayInfo.resonances.size()];
        calculators[i]  = new SpecialResonanceCalculator(parameters, i);
        integrals[i]    = new fpcomplex *[decayInfo.resonances.size()];

        for(int j = 0; j < decayInfo.resonances.size(); ++j) {
            integrals[i][j]   = new fpcomplex(0, 0);
            integrators[i][j] = new SpecialResonanceIntegrator(parameters, i, j);
        }
    }

    setSeparateNorm();
}

void Amp3Body::populateArrays() {
    PdfBase::populateArrays();

    // save our efficiency function.  Resonance's are saved first, then the efficiency function.  Take -1 as efficiency!
    efficiencyFunction = host_function_table.size() - 1;
}
__host__ void Amp3Body::setDataSize(unsigned int dataSize, unsigned int evtSize) {
    // Default 3 is m12, m13, evtNum
    totalEventSize = evtSize;
    if(totalEventSize < 3)
        throw GooFit::GeneralError("totalEventSize {} must be 3 or more", totalEventSize);

    // if (cachedWaves) delete cachedWaves;
    if(cachedWaves[0]) {
        for(auto &cachedWave : cachedWaves) {
            delete cachedWave;
            cachedWave = nullptr;
        }
    }

    numEntries = dataSize;

    for(int i = 0; i < 16; i++) {
#ifdef GOOFIT_MPI
        cachedWaves[i] = new thrust::device_vector<fpcomplex>(m_iEventsPerTask);
#else
        cachedWaves[i] = new thrust::device_vector<fpcomplex>(dataSize);
#endif
        void *dummy = thrust::raw_pointer_cast(cachedWaves[i]->data());
        MEMCPY_TO_SYMBOL(cResonances, &dummy, sizeof(fpcomplex *), i * sizeof(fpcomplex *), cudaMemcpyHostToDevice);
    }

    setForceIntegrals();
}

__host__ fptype Amp3Body::normalize() {
    recursiveSetNormalization(1.0); // Not going to normalize efficiency,
    // so set normalization factor to 1 so it doesn't get multiplied by zero.
    // Copy at this time to ensure that the SpecialResonanceCalculators, which need the efficiency,
    // don't get zeroes through multiplying by the normFactor.
    // we need to update the normal here, as values are used at this point.

    host_normalizations.sync(d_normalizations);

    int totalBins = _m12.getNumBins() * _m13.getNumBins();

    if(!dalitzNormRange) {
        gooMalloc((void **)&dalitzNormRange, 6 * sizeof(fptype));
    }

    // This line runs once
    static std::array<fptype, 6> host_norms{{0, 0, 0, 0, 0, 0}};

    std::array<fptype, 6> current_host_norms{{_m12.getLowerLimit(),
                                              _m12.getUpperLimit(),
                                              static_cast<fptype>(_m12.getNumBins()),
                                              _m13.getLowerLimit(),
                                              _m13.getUpperLimit(),
                                              static_cast<fptype>(_m13.getNumBins())}};

    if(host_norms != current_host_norms) {
        host_norms = current_host_norms;
        MEMCPY(dalitzNormRange, host_norms.data(), 6 * sizeof(fptype), cudaMemcpyHostToDevice);
    }

    for(unsigned int i = 0; i < decayInfo.resonances.size(); ++i) {
        redoIntegral[i] = forceRedoIntegrals;

        if(!(decayInfo.resonances[i]->parametersChanged()))
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

    for(int i = 0; i < decayInfo.resonances.size(); ++i) {
        // grab the index for this resonance.
        calculators[i]->setResonanceIndex(decayInfo.resonances[i]->getFunctionIndex());
        calculators[i]->setDalitzIndex(getFunctionIndex());
        if(redoIntegral[i]) {
#ifdef GOOFIT_MPI
            thrust::transform(
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex + m_iEventsPerTask, arrayAddress, eventSize)),
                strided_range<thrust::device_vector<fpcomplex>::iterator>(
                    cachedWaves[i]->begin(), cachedWaves[i]->end(), 1)
                    .begin(),
                *(calculators[i]));
#else
            thrust::transform(
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, arrayAddress, eventSize)),
                strided_range<thrust::device_vector<fpcomplex>::iterator>(
                    cachedWaves[i]->begin(), cachedWaves[i]->end(), 1)
                    .begin(),
                *(calculators[i]));
#endif
        }

        // Possibly this can be done more efficiently by exploiting symmetry?
        for(int j = 0; j < decayInfo.resonances.size(); ++j) {
            if((!redoIntegral[i]) && (!redoIntegral[j]))
                continue;

            integrators[i][j]->setDalitzIndex(getFunctionIndex());
            integrators[i][j]->setResonanceIndex(decayInfo.resonances[i]->getFunctionIndex());
            // integrators[i][j]->setEfficiencyIndex(efficiencyFunction);
            integrators[i][j]->setEfficiencyIndex(decayInfo.resonances[j]->getFunctionIndex());
            thrust::constant_iterator<int> effFunc(efficiencyFunction);
            fpcomplex dummy(0, 0);
            thrust::plus<fpcomplex> complexSum;
            (*(integrals[i][j])) = thrust::transform_reduce(
                thrust::make_zip_iterator(thrust::make_tuple(binIndex, arrayAddress, effFunc)),
                thrust::make_zip_iterator(thrust::make_tuple(binIndex + totalBins, arrayAddress, effFunc)),
                *(integrators[i][j]),
                dummy,
                complexSum);
        }
    }

    // End of time-consuming integrals.
    fpcomplex sumIntegral(0, 0);

    for(unsigned int i = 0; i < decayInfo.resonances.size(); ++i) {
        // int param_i = parameters + resonanceOffset_DP + resonanceSize * i;
        fpcomplex amplitude_i(host_parameters[parametersIdx + i * 2 + 1], host_parameters[parametersIdx + i * 2 + 2]);

        // printf("i:%i - %f,%f\n", i, amplitude_i.real(), amplitude_i.imag());

        for(unsigned int j = 0; j < decayInfo.resonances.size(); ++j) {
            // int param_j = parameters + resonanceOffset_DP + resonanceSize * j;
            fpcomplex amplitude_j(host_parameters[parametersIdx + j * 2 + 1],
                                  -host_parameters[parametersIdx + j * 2 + 2]);

            // printf("j:%i - %f,%f\n", j, amplitude_j.real(), amplitude_j.imag());
            // Notice complex conjugation
            // printf("%f %f %f %f %f %f\n", amplitude_i.real(), amplitude_i.imag(), amplitude_j.real(),
            // amplitude_j.imag(), (*(integrals[i][j])).real(), (*(integrals[i][j])).imag() );
            sumIntegral += amplitude_i * amplitude_j * (*(integrals[i][j]));
        }
    }

    fptype ret           = sumIntegral.real(); // That complex number is a square, so it's fully real
    double binSizeFactor = 1;
    binSizeFactor *= _m12.getBinSize();
    binSizeFactor *= _m13.getBinSize();
    ret *= binSizeFactor;

    host_normalizations[normalIdx + 1] = 1.0 / ret;
    cachedNormalization                = 1.0 / ret;
    // printf("%f %f\n", ret, binSizeFactor);
    return ret;
}

__host__ fpcomplex Amp3Body::sumCachedWave(size_t i) const {
    const thrust::device_vector<fpcomplex> &vec = getCachedWaveNoCopy(i);

    fpcomplex ret = thrust::reduce(vec.begin(), vec.end(), fpcomplex(0, 0), thrust::plus<fpcomplex>());

    return ret;
}

__host__ const std::vector<std::complex<fptype>> Amp3Body::getCachedWave(size_t i) const {
    // TODO: This calls itself immediatly ?
    auto ret_thrust = getCachedWave(i);
    std::vector<std::complex<fptype>> ret(ret_thrust.size());
    thrust::copy(ret_thrust.begin(), ret_thrust.end(), ret.begin());
    return ret;
}

__host__ std::vector<std::vector<fptype>> Amp3Body::fit_fractions() {
    GOOFIT_DEBUG("Performing fit fraction calculation, should already have a cache (does not use normalization grid)");

    size_t n_res    = getDecayInfo().resonances.size();
    size_t nEntries = getCachedWaveNoCopy(0).size();

    std::vector<fpcomplex> coefs(n_res);
    std::transform(getDecayInfo().resonances.begin(),
                   getDecayInfo().resonances.end(),
                   coefs.begin(),
                   [](ResonancePdf *res) { return fpcomplex(res->amp_real.getValue(), res->amp_imag.getValue()); });

    fptype buffer_all = 0;
    fptype buffer;
    fpcomplex coef_i;
    fpcomplex coef_j;
    fpcomplex cached_i_val;
    fpcomplex cached_j_val;

    thrust::device_vector<fpcomplex> cached_i;
    thrust::device_vector<fpcomplex> cached_j;
    std::vector<std::vector<fptype>> Amps_int(n_res, std::vector<fptype>(n_res));

    for(size_t i = 0; i < n_res; i++) {
        for(size_t j = 0; j < n_res; j++) {
            buffer   = 0;
            cached_i = getCachedWaveNoCopy(i);
            cached_j = getCachedWaveNoCopy(j);
            coef_i   = coefs[i];
            coef_j   = coefs[j];

            buffer += thrust::transform_reduce(
                thrust::make_zip_iterator(thrust::make_tuple(cached_i.begin(), cached_j.begin())),
                thrust::make_zip_iterator(thrust::make_tuple(cached_i.end(), cached_j.end())),
                CoefSumFunctor(coef_i, coef_j),
                (fptype)0.0,
                thrust::plus<fptype>());

            buffer_all += buffer;
            Amps_int[i][j] = (buffer / nEntries);
        }
    }

    fptype total_PDF = buffer_all / nEntries;

    std::vector<std::vector<fptype>> ff(n_res, std::vector<fptype>(n_res));

    for(size_t i = 0; i < n_res; i++)
        for(size_t j = 0; j < n_res; j++)
            ff[i][j] = (Amps_int[i][j] / total_PDF);

    return ff;
}

} // namespace GooFit
