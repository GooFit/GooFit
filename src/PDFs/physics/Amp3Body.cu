#include <mcbooster/Evaluate.h>
#include <mcbooster/EvaluateArray.h>
#include <mcbooster/GContainers.h>
#include <mcbooster/GFunctional.h>
#include <mcbooster/GTypes.h>
#include <mcbooster/Generate.h>
#include <mcbooster/Vector4R.h>

#include <goofit/Error.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/Amp3Body.h>
#include <goofit/PDFs/physics/Amp3BodyBase.h>
#include <goofit/PDFs/physics/detail/Dim2.h>
#include <goofit/PDFs/physics/detail/SpecialResonanceCalculator.h>
#include <goofit/PDFs/physics/detail/SpecialResonanceIntegrator.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>
#include <goofit/detail/Complex.h>

#include <thrust/copy.h>
#include <thrust/transform_reduce.h>

#include <array>
#include <vector>

namespace GooFit {

// Functor used for fit fraction sum
struct CoefSumFunctor {
    fpcomplex coef_i;
    fpcomplex coef_j;

    CoefSumFunctor(fpcomplex coef_i, fpcomplex coef_j)
        : coef_i(coef_i)
        , coef_j(coef_j) {}

    __device__ auto operator()(thrust::tuple<fpcomplex, fpcomplex> val) -> fptype {
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
// this needs to be large enough to hold all samples
__device__ fpcomplex *cResonances[16 * 20];

__device__ inline auto parIndexFromResIndex_DP(int resIndex) -> int {
    return resonanceOffset_DP + resIndex * resonanceSize;
}

__device__ auto device_DalitzPlot(fptype *evt, ParameterContainer &pc) -> fptype {
    int num_obs = pc.getNumObservables();
    int id_m12  = pc.getObservable(0);
    int id_m13  = pc.getObservable(1);
    int id_num  = pc.getObservable(2);

    fptype m12 = RO_CACHE(evt[id_m12]);
    fptype m13 = RO_CACHE(evt[id_m13]);

    unsigned int numResonances = pc.getConstant(0);
    unsigned int cacheToUse    = pc.getConstant(1);

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
        fpcomplex me = RO_CACHE(cResonances[i + (16 * cacheToUse)][evtNum]);

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

int Amp3Body::cacheCount                         = 0;
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
    MEMCPY_TO_SYMBOL(c_mother_meson_radius, &decay.mother_meson_radius, sizeof(fptype), 0, cudaMemcpyHostToDevice);

    // registered to 0 position
    registerConstant(decayInfo.resonances.size());

    cacheToUse = cacheCount++;
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
__host__ void Amp3Body::setDataSize(unsigned int dataSize, unsigned int evtSize, unsigned int offset) {
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

    numEntries  = dataSize;
    eventOffset = offset;

    for(int i = 0; i < 16; i++) {
#ifdef GOOFIT_MPI
        cachedWaves[i] = new thrust::device_vector<fpcomplex>(m_iEventsPerTask);
#else
        cachedWaves[i] = new thrust::device_vector<fpcomplex>(dataSize);
#endif
        void *dummy = thrust::raw_pointer_cast(cachedWaves[i]->data());
        MEMCPY_TO_SYMBOL(cResonances,
                         &dummy,
                         sizeof(fpcomplex *),
                         ((16 * cacheToUse) + i) * sizeof(fpcomplex *),
                         cudaMemcpyHostToDevice);
    }

    setForceIntegrals();
}

__host__ auto Amp3Body::normalize() -> fptype {
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
    }
    MEMCPY(dalitzNormRange, host_norms.data(), 6 * sizeof(fptype), cudaMemcpyHostToDevice);
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
    thrust::counting_iterator<int> eventIndex(eventOffset);

    for(int i = 0; i < decayInfo.resonances.size(); ++i) {
        // grab the index for this resonance.
        calculators[i]->setResonanceIndex(decayInfo.resonances[i]->getFunctionIndex());
        calculators[i]->setDalitzIndex(getFunctionIndex());
        if(redoIntegral[i]) {
#ifdef GOOFIT_MPI
            thrust::transform(
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex + m_iEventsPerTask, dataArray, eventSize)),
                strided_range<thrust::device_vector<fpcomplex>::iterator>(
                    cachedWaves[i]->begin(), cachedWaves[i]->end(), 1)
                    .begin(),
                *(calculators[i]));
#else
            thrust::transform(
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                // was this correct before?
                // thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, dataArray, eventSize)),
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

        for(unsigned int j = 0; j < decayInfo.resonances.size(); ++j) {
            // int param_j = parameters + resonanceOffset_DP + resonanceSize * j;
            fpcomplex amplitude_j(host_parameters[parametersIdx + j * 2 + 1],
                                  -host_parameters[parametersIdx + j * 2 + 2]);

            // Notice complex conjugation
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
    return ret;
}

__host__ auto Amp3Body::sumCachedWave(size_t i) const -> fpcomplex {
    const thrust::device_vector<fpcomplex> &vec = getCachedWaveNoCopy(i);

    fpcomplex ret = thrust::reduce(vec.begin(), vec.end(), fpcomplex(0, 0), thrust::plus<fpcomplex>());

    return ret;
}

__host__ auto Amp3Body::getCachedWave(size_t i) const -> const std::vector<std::complex<fptype>> {
    // TODO: This calls itself immediately ?
    auto ret_thrust = getCachedWave(i);
    std::vector<std::complex<fptype>> ret(ret_thrust.size());
    thrust::copy(ret_thrust.begin(), ret_thrust.end(), ret.begin());
    return ret;
}

__host__ auto Amp3Body::fit_fractions() -> std::vector<std::vector<fptype>> {
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

__host__ auto Amp3Body::GenerateSig(unsigned int numEvents, int seed) -> std::
    tuple<mcbooster::ParticlesSet_h, mcbooster::VariableSet_h, mcbooster::RealVector_h, mcbooster::RealVector_h> {
    // Must configure our functions before any calculations!
    // setupObservables();
    // setIndices();

    initialize();

    // Defining phase space
    std::vector<mcbooster::GReal_t> masses{decayInfo.daug1Mass, decayInfo.daug2Mass, decayInfo.daug3Mass};
    mcbooster::PhaseSpace phsp(decayInfo.motherMass, masses, numEvents, generation_offset);

    if(seed != 0) {
        phsp.SetSeed(seed);
    } else {
        GOOFIT_INFO("Current generator seed {}, offset {}", phsp.GetSeed(), generation_offset);
    }

    // Generating numEvents events. Events are all generated inside the phase space with uniform distribution in
    // momentum space. Events must be weighted to have phase space distribution
    phsp.Generate(mcbooster::Vector4R(decayInfo.motherMass, 0.0, 0.0, 0.0));

    auto d1 = phsp.GetDaughters(0);
    auto d2 = phsp.GetDaughters(1);
    auto d3 = phsp.GetDaughters(2);

    mcbooster::ParticlesSet_d pset(3);
    pset[0] = &d1;
    pset[1] = &d2;
    pset[2] = &d3;

    auto SigGen_M12_d = mcbooster::RealVector_d(numEvents);
    auto SigGen_M13_d = mcbooster::RealVector_d(numEvents);
    auto SigGen_M23_d = mcbooster::RealVector_d(numEvents);

    mcbooster::VariableSet_d VarSet_d(3);
    VarSet_d[0] = &SigGen_M12_d;
    VarSet_d[1] = &SigGen_M23_d;
    VarSet_d[2] = &SigGen_M13_d;

    // Evaluating invariant masses for each event
    Dim2 eval = Dim2();
    mcbooster::EvaluateArray<Dim2>(eval, pset, VarSet_d);

    mcbooster::VariableSet_d GooVarSet_d(3);
    GooVarSet_d[0] = VarSet_d[0];
    GooVarSet_d[1] = VarSet_d[2];
    GooVarSet_d[2] = VarSet_d[1];

    auto h1 = new mcbooster::Particles_h(d1);
    auto h2 = new mcbooster::Particles_h(d2);
    auto h3 = new mcbooster::Particles_h(d3);

    mcbooster::ParticlesSet_h ParSet(3);
    ParSet[0] = h1;
    ParSet[1] = h2;
    ParSet[2] = h3;

    auto SigGen_M12_h = new mcbooster::RealVector_h(SigGen_M12_d);
    auto SigGen_M23_h = new mcbooster::RealVector_h(SigGen_M23_d);
    auto SigGen_M13_h = new mcbooster::RealVector_h(SigGen_M13_d);

    mcbooster::VariableSet_h VarSet(3);
    VarSet[0] = SigGen_M12_h;
    VarSet[1] = SigGen_M23_h;
    VarSet[2] = SigGen_M13_h;

    mcbooster::RealVector_d weights(phsp.GetWeights());
    phsp.FreeResources();

    auto DS = new mcbooster::RealVector_d(3 * numEvents);
    thrust::counting_iterator<int> eventNumber(0);

#pragma unroll

    for(int i = 0; i < 2; ++i) {
        mcbooster::strided_range<mcbooster::RealVector_d::iterator> sr(DS->begin() + i, DS->end(), 3);
        thrust::copy(GooVarSet_d[i]->begin(), GooVarSet_d[i]->end(), sr.begin());
    }

    mcbooster::strided_range<mcbooster::RealVector_d::iterator> sr(DS->begin() + 2, DS->end(), 3);
    thrust::copy(eventNumber, eventNumber + numEvents, sr.begin());

    // Giving events to GooFit. Format of dev_evt_array must be (s12, s13, eventNumber). s23 is calculated automatically
    // in src/PDFs/physics/detail/SpecialResonanceCalculator.cu
    dev_event_array = thrust::raw_pointer_cast(DS->data());
    setDataSize(numEvents, 3);

    generation_no_norm = true; // we need no normalization for generation, but we do need to make sure that norm = 1;
    SigGenSetIndices();
    copyParams();
    normalize();
    setForceIntegrals();
    host_normalizations.sync(d_normalizations);

    auto fc = fitControl;
    setFitControl(std::make_shared<ProbFit>());

    thrust::device_vector<fptype> results;
    GooPdf::evaluate_with_metric(results);

    // evaluating amplitudes for generated events, amplitudes are incorporated in weights
    thrust::transform(
        results.begin(), results.end(), weights.begin(), weights.begin(), thrust::multiplies<mcbooster::GReal_t>());

    // Filing accept/reject flags for resonant distribution for each generated event
    mcbooster::BoolVector_d flags(numEvents);
    fillMCFlags(flags, weights, numEvents);

    auto weights_h = mcbooster::RealVector_h(weights);
    auto results_h = mcbooster::RealVector_h(results);
    auto flags_h   = mcbooster::BoolVector_h(flags);
    cudaDeviceSynchronize();

    setFitControl(fc);

    return std::make_tuple(ParSet, VarSet, weights_h, flags_h);
}
} // namespace GooFit
