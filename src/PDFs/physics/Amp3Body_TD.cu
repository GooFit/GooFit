#include <goofit/Error.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/Amp3Body_TD.h>
#include <goofit/PDFs/physics/detail/SpecialComplexSum.h>
#include <goofit/PDFs/physics/detail/SpecialDalitzIntegrator.h>
#include <goofit/PDFs/physics/detail/SpecialWaveCalculator.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

#include <thrust/transform_reduce.h>

#ifdef GOOFIT_MPI
#include <mpi.h>
#endif

namespace GooFit {

const int resonanceOffset = 10; // Offset of the first resonance into the parameter index array
// Offset is number of parameters, constant index, indices for tau, xmix, and ymix, index
// of resolution function, and finally number of resonances (not calculable from nP
// because we don't know what the efficiency and time resolution might need). Efficiency
// and time-resolution parameters are after the resonance information.

const unsigned int SPECIAL_RESOLUTION_FLAG = 999999999;

// The function of this array is to hold all the cached waves; specific
// waves are recalculated when the corresponding resonance mass or width
// changes. Note that in a multithread environment each thread needs its
// own cache, hence the '10'. Ten threads should be enough for anyone!

// NOTE: only one set of wave holders is supported currently.
// this needs to be large enough to hold all samples
__device__ WaveHolder_s *cWaves[16 * 20];

__device__ inline auto parIndexFromResIndex(int resIndex) -> int { return resonanceOffset + resIndex * resonanceSize; }

__device__ auto getResonanceAmplitude(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) -> fpcomplex {
    auto func = reinterpret_cast<resonance_function_ptr>(d_function_table[pc.funcIdx]);
    return (*func)(m12, m13, m23, pc);
}

__device__ auto device_Tddp(fptype *evt, ParameterContainer &pc) -> fptype {
    int num_parameters  = pc.getNumParameters();
    int num_constants   = pc.getNumConstants();
    int num_observables = pc.getNumObservables();

    int id_m12 = pc.getObservable(2);
    int id_m13 = pc.getObservable(3);
    int id_num = pc.getObservable(4);
    int id_mis = 0;
    int id_tag = 0;
    if(num_observables > 5) {
        id_mis = pc.getObservable(5);
    }

    fptype m12 = RO_CACHE(evt[id_m12]);
    fptype m13 = RO_CACHE(evt[id_m13]);

    unsigned int numResonances = pc.getConstant(0);
    // int numResonances = 1;

    if(!inDalitz(m12, m13, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass)) {
        unsigned int endEfficiencyFunc = pc.getConstant(3);
        pc.incrementIndex(1, num_parameters, num_constants, num_observables, 1);

        // increment the resonances
        for(int i = 0; i < numResonances; i++)
            pc.incrementIndex();

        // increment the resolution function
        pc.incrementIndex();

        // increment our efficiency function
        // pc.incrementIndex();
        while(pc.funcIdx < endEfficiencyFunc)
            pc.incrementIndex();
        return 0;
    }

    auto evtNum = static_cast<int>(floor(0.5 + RO_CACHE(evt[id_num])));

    fpcomplex sumWavesA(0, 0);
    fpcomplex sumWavesB(0, 0);
    fpcomplex sumRateAA(0, 0);
    fpcomplex sumRateAB(0, 0);
    fpcomplex sumRateBB(0, 0);

    unsigned int cacheToUse = pc.getConstant(1);
    fptype mistag           = 0;
    // fptype mistag = pc.getConstant(2);

    for(int i = 0; i < numResonances; ++i) {
        // int paramIndex = parIndexFromResIndex(i);
        fpcomplex amp{pc.getParameter(i * 2 + 5), pc.getParameter(i * 2 + 6)};

        // fpcomplex matrixelement(thrust::get<0>(cWaves[cacheToUse][evtNum*numResonances + i]),
        //				     thrust::get<1>(cWaves[cacheToUse][evtNum*numResonances + i]));
        // Note, to make this more efficient we should change it to only an array of fptype's, and read double2 at a
        // time.
        int index_cWave = i + (16 * cacheToUse);
        fpcomplex ai{RO_CACHE(cWaves[index_cWave][evtNum].ai_real), RO_CACHE(cWaves[index_cWave][evtNum].ai_imag)};
        fpcomplex bi{RO_CACHE(cWaves[index_cWave][evtNum].bi_real), RO_CACHE(cWaves[index_cWave][evtNum].bi_imag)};

        fpcomplex matrixelement = ai * amp;
        sumWavesA += matrixelement;

        // matrixelement = fpcomplex(thrust::get<2>(cWaves[cacheToUse][evtNum*numResonances + i]),
        //				       thrust::get<3>(cWaves[cacheToUse][evtNum*numResonances + i]));
        matrixelement = bi * amp;
        sumWavesB += matrixelement;
    }

    int id_time  = pc.getObservable(0);
    int id_sigma = pc.getObservable(1);
    if(num_observables > 6) {
        id_tag = pc.getObservable(6);
    }
    fptype _tau = pc.getParameter(0);
    // fptype _xmixing = pc.getParameter(1);
    // fptype _ymixing = pc.getParameter(2);
    fptype _xmixing0 = pc.getParameter(1);
    fptype _ymixing0 = pc.getParameter(2);
    fptype _deltax   = pc.getParameter(3);
    fptype _deltay   = pc.getParameter(4);
    fptype _xmixing  = 0;
    fptype _ymixing  = 0;
    // int _charmtag = evt[id_tag];
    // auto _charmtag = static_cast<int>(floor(0.5 + RO_CACHE(evt[id_tag])));
    int _charmtag = 0;
    if(num_observables > 6) {
        _charmtag = RO_CACHE(evt[id_tag]);
    }
    if(_charmtag == 1) {
        _xmixing = _xmixing0 + _deltax;
        _ymixing = _ymixing0 + _deltay;
    } else if(_charmtag == -1) {
        _xmixing = _xmixing0 - _deltax;
        _ymixing = _ymixing0 - _deltay;
    } else {
        _xmixing = _xmixing0;
        _ymixing = _ymixing0;
    }

    fptype _time  = RO_CACHE(evt[id_time]);
    fptype _sigma = RO_CACHE(evt[id_sigma]);

    // TODO: Test that we have a special flag by comparing size of numconstants?
    // fptype special_flag = pc.getConstant(3);

    // if ((gpuDebug & 1) && (0 == BLOCKIDX) && (0 == THREADIDX))
    // if (0 == evtNum) printf("TDDP: (%f, %f) (%f, %f)\n", sumWavesA.real, sumWavesA.imag, sumWavesB.real,
    // sumWavesB.imag);
    // printf("TDDP: %f %f %f %f | %f %f %i\n", m12, m13, _time, _sigma, _xmixing, _tau, evtNum);

    /*
    fptype ret = 0;
    ret += (norm2(sumWavesA) + norm2(sumWavesB))*cosh(_ymixing * _time);
    ret += (norm2(sumWavesA) - norm2(sumWavesB))*cos (_xmixing * _time);
    sumWavesA *= conj(sumWavesB);
    ret -= 2*sumWavesA.real * sinh(_ymixing * _time);
    ret -= 2*sumWavesA.imag * sin (_xmixing * _time); // Notice sign difference wrt to Mikhail's code, because I have
    AB* and he has A*B.
    ret *= exp(-_time);
    */

    fptype term1 = thrust::norm(sumWavesA) + thrust::norm(sumWavesB);
    fptype term2 = thrust::norm(sumWavesA) - thrust::norm(sumWavesB);
    sumWavesA *= conj(sumWavesB);
    // printf("(%i, %i) TDDP: %f %f %f %f %f %f %f\n", BLOCKIDX, THREADIDX, term1, term2, sumWavesA.real,
    // sumWavesA.imag, m12, m13, _tau);

    // Cannot use callFunction on resolution function.
    // int effFunctionIdx = parIndexFromResIndex(numResonances);
    // int resFunctionIdx = RO_CACHE(indices[5]);
    // int resFunctionPar = 2 + effFunctionIdx;
    fptype ret = 0;
    // int md0_offset     = 0;

    // if(resFunctionIdx == SPECIAL_RESOLUTION_FLAG) {
    // In this case there are multiple resolution functions, they are stored after the efficiency function,
    // and which one we use depends on the measured mother-particle mass.
    //    md0_offset     = 1;
    // int id_massd0 = pc.constants[pc.constantIdx + 6];
    // fptype massd0  = RO_CACHE(evt[id_massd0]);
    // fptype minMass = RO_CACHE(pc.constants[pc.constantIdx + 7]);
    // fptype md0Step = RO_CACHE(pc.constants[pc.constantIdx + 8]);
    // int res_to_use = (massd0 <= minMass) ? 0 : static_cast<int>(floor((massd0 - minMass) / md0Step));
    // int maxFcn     = RO_CACHE(indices[2 + effFunctionIdx]);

    //    if(res_to_use > maxFcn)
    //        res_to_use = maxFcn;

    // Now calculate index of resolution function.
    // At the end of the array are indices efficiency_function, efficiency_parameters, maxFcn, res_function_1,
    // res_function_1_nP, par1, par2 ... res_function_2, res_function_2_nP, ...
    // res_to_use = 3 + effFunctionIdx + res_to_use * (2 + RO_CACHE(indices[effFunctionIdx + 4]));
    // NB this assumes all resolution functions have the same number of parameters. The number
    // of parameters in the first resolution function is stored in effFunctionIdx+3; add one to
    // account for the index of the resolution function itself in the device function table, one
    // to account for the number-of-parameters index, and this is the space taken up by each
    // resolution function. Multiply by res_to_use to get the number of spaces to skip to get to
    // the one we want.

    // resFunctionIdx = RO_CACHE(indices[res_to_use]);
    // resFunctionPar = res_to_use + 1;
    //}

    pc.incrementIndex(1, num_parameters, num_constants, num_observables, 1);

    // increment over resonance functions here?
    for(int i = 0; i < numResonances; i++)
        pc.incrementIndex();

    ret = (*(reinterpret_cast<device_resfunction_ptr>(d_function_table[pc.funcIdx])))(
        term1, term2, sumWavesA.real(), sumWavesA.imag(), _tau, _time, _xmixing, _ymixing, _sigma, pc);

    // For the reversed (mistagged) fraction, we make the
    // interchange A <-> B. So term1 stays the same,
    // term2 changes sign, and AB* becomes BA*.
    // Efficiency remains the same for the mistagged part,
    // because it depends on the momenta of the pi+ and pi-,
    // which don't change even though we tagged a D0 as D0bar.

    // fptype mistag = RO_CACHE(functorConstants[RO_CACHE(indices[1]) + 5]);
    mistag = 0;
    if(num_observables > 5) {
        mistag = evt[id_mis];
    } // mistag = 0;

    fptype xfix   = 0.0039;
    fptype yfix   = 0.0065;
    fptype taufix = 0.4101;

    if(mistag > 0) { // This should be either true or false for all events, so no branch is caused.
        // See header file for explanation of 'mistag' variable - it is actually the probability
        // of having the correct sign, given that we have a correctly reconstructed D meson.
        // mistag = evt[id_mis];
        ret *= (1 - 2 * mistag);
        // The following formats differently in clang-format 8
        // clang-format off
        ret += mistag
               * (*(reinterpret_cast<device_resfunction_ptr>(d_function_table[pc.funcIdx])))(
                    term1, -term2, sumWavesA.real(), -sumWavesA.imag(), taufix, _time, xfix, yfix, _sigma, pc);
        ret += mistag
               * (*(reinterpret_cast<device_resfunction_ptr>(d_function_table[pc.funcIdx])))(
                    term1, term2, sumWavesA.real(), sumWavesA.imag(), taufix, _time, xfix, yfix, _sigma, pc);

        // clang-format on
    }

    // increment our resolution function
    pc.incrementIndex();

    fptype eff = callFunction(evt, pc);
    // internalDebug = 0;
    ret *= eff;

    return ret;
}

int Amp3Body_TD::cacheCount                = 0;
__device__ device_function_ptr ptr_to_Tddp = device_Tddp;

__host__ Amp3Body_TD::Amp3Body_TD(std::string n,
                                  Observable _dtime,
                                  Observable _sigmat,
                                  Observable m12,
                                  Observable m13,
                                  EventNumber eventNumber,
                                  DecayInfo3t decay,
                                  MixingTimeResolution *r,
                                  GooPdf *efficiency,
                                  Observable *mistag,
                                  Observable *charmtag)
    : Amp3BodyBase("Amp3Body_TD", n, _dtime, _sigmat, m12, m13, eventNumber)
    , decayInfo(decay)
    , _m12(m12)
    , _m13(m13)
    , resolution(r)
    , _mistag(*mistag)
    , _charmtag(*charmtag)
    , totalEventSize(6) // Default 5 = m12, m13, time, sigma_t, evtNum
{
    for(auto &cachedWave : cachedWaves)
        cachedWave = nullptr;

    registerObservable(*mistag);
    totalEventSize++;

    registerObservable(*charmtag);
    totalEventSize++;

    MEMCPY_TO_SYMBOL(c_motherMass, &decay.motherMass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug1Mass, &decay.daug1Mass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug2Mass, &decay.daug2Mass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug3Mass, &decay.daug3Mass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_meson_radius, &decay.meson_radius, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_mother_meson_radius, &decay.mother_meson_radius, sizeof(fptype), 0, cudaMemcpyHostToDevice);

    registerParameter(decay._tau);
    registerParameter(decay._xmixing);
    registerParameter(decay._ymixing);
    registerParameter(decay._deltax);
    registerParameter(decay._deltay);

    setD0Fraction(0.5);

    if(resolution->getDeviceFunction() < 0)
        throw GooFit::GeneralError("The resolution device function index {} must be more than 0",
                                   resolution->getDeviceFunction());

    registerConstant(decay.resonances.size());

    cacheToUse = cacheCount++;
    registerConstant(cacheToUse);

    if(mistag == nullptr)
        registerConstant(1);
    else
        registerConstant(0);

    for(auto &resonance : decay.resonances) {
        registerParameter(resonance->amp_real);
        registerParameter(resonance->amp_imag);
        components.push_back(resonance);
    }

    components.push_back(resolution);
    components.push_back(efficiency);

    // this is the funcID after the efficiency routine
    registerConstant(0);

    registerFunction("ptr_to_Tddp", ptr_to_Tddp);

    initialize();

    redoIntegral = new bool[decay.resonances.size()];
    cachedMasses = new fptype[decay.resonances.size()];
    cachedWidths = new fptype[decay.resonances.size()];
    integrals    = new ThreeComplex **[decay.resonances.size()];
    integrators  = new SpecialDalitzIntegrator **[decay.resonances.size()];
    calculators  = new SpecialWaveCalculator *[decay.resonances.size()];

    for(int i = 0; i < decay.resonances.size(); ++i) {
        redoIntegral[i] = true;
        cachedMasses[i] = -1;
        cachedWidths[i] = -1;
        integrators[i]  = new SpecialDalitzIntegrator *[decay.resonances.size()];
        calculators[i]  = new SpecialWaveCalculator(parameters, i);
        integrals[i]    = new ThreeComplex *[decay.resonances.size()];

        for(int j = 0; j < decay.resonances.size(); ++j) {
            integrals[i][j]   = new ThreeComplex(0, 0, 0, 0, 0, 0);
            integrators[i][j] = new SpecialDalitzIntegrator(parameters, i, j);
        }
    }

    setSeparateNorm();
}

__host__ Amp3Body_TD::Amp3Body_TD(std::string n,
                                  Observable _dtime,
                                  Observable _sigmat,
                                  Observable m12,
                                  Observable m13,
                                  EventNumber eventNumber,
                                  DecayInfo3t decay,
                                  std::vector<MixingTimeResolution *> &r,
                                  GooPdf *efficiency,
                                  Observable md0,
                                  Observable *mistag,
                                  Observable *charmtag)
    : Amp3BodyBase("Amp3Body_TD", n, _dtime, _sigmat, m12, m13, eventNumber, md0)
    , decayInfo(decay)
    , _m12(m12)
    , _m13(m13)
    , resolution(
          r[0]) // Only used for normalization, which only depends on x and y - it doesn't matter which one we use.
    , _mistag(*mistag)
    , _charmtag(*charmtag)
    , totalEventSize(6) // This case adds the D0 mass by default.
{
    for(auto &cachedWave : cachedWaves)
        cachedWave = nullptr;

    registerObservable(*mistag);
    totalEventSize++;

    registerObservable(*charmtag);
    totalEventSize++;

    MEMCPY_TO_SYMBOL(c_motherMass, &decay.motherMass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug1Mass, &decay.daug1Mass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug2Mass, &decay.daug2Mass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug3Mass, &decay.daug3Mass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_meson_radius, &decay.meson_radius, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_mother_meson_radius, &decay.mother_meson_radius, sizeof(fptype), 0, cudaMemcpyHostToDevice);

    registerParameter(decay._tau);
    registerParameter(decay._xmixing);
    registerParameter(decay._ymixing);
    registerParameter(decay._deltax);
    registerParameter(decay._deltay);
    printf("Multiple resolution functions not supported yet!\n");

    registerConstant(decayInfo.resonances.size());

    cacheToUse = cacheCount++;
    registerConstant(cacheToUse);

    if(mistag)
        registerConstant(1);
    else
        registerConstant(0);

    registerConstant(SPECIAL_RESOLUTION_FLAG);

    // TODO: Do these need to be set as constants?
    registerConstant(md0.getLowerLimit());
    registerConstant((md0.getUpperLimit() - md0.getLowerLimit()) / r.size());

    for(auto &resonance : decayInfo.resonances) {
        registerParameter(resonance->amp_real);
        registerParameter(resonance->amp_imag);

        components.push_back(resonance);
    }

    // components.push_back(resolution);

    for(auto &i : r) {
        if(i->getDeviceFunction() < 0)
            throw GooFit::GeneralError("Device function index {} must be more than 0", i->getDeviceFunction());

        components.push_back(i);
    }

    components.push_back(efficiency);

    registerFunction("ptr_to_Tddp", ptr_to_Tddp);

    initialize();

    // this is the funcID after the efficiency routine
    registerConstant(0);

    redoIntegral = new bool[decay.resonances.size()];
    cachedMasses = new fptype[decay.resonances.size()];
    cachedWidths = new fptype[decay.resonances.size()];
    integrals    = new ThreeComplex **[decay.resonances.size()];
    integrators  = new SpecialDalitzIntegrator **[decay.resonances.size()];
    calculators  = new SpecialWaveCalculator *[decay.resonances.size()];

    for(int i = 0; i < decay.resonances.size(); ++i) {
        redoIntegral[i] = true;
        cachedMasses[i] = -1;
        cachedWidths[i] = -1;
        integrators[i]  = new SpecialDalitzIntegrator *[decay.resonances.size()];
        calculators[i]  = new SpecialWaveCalculator(parameters, i);
        integrals[i]    = new ThreeComplex *[decay.resonances.size()];

        for(int j = 0; j < decay.resonances.size(); ++j) {
            integrals[i][j]   = new ThreeComplex(0, 0, 0, 0, 0, 0);
            integrators[i][j] = new SpecialDalitzIntegrator(parameters, i, j);
        }
    }

    setSeparateNorm();
}

// Note: We need to manually populate the arrays so we can track the efficiency function!
__host__ void Amp3Body_TD::populateArrays() {
    // populate all the arrays
    GOOFIT_TRACE("Amp3Body_TD: Populating Arrays for {}", getName());

    // reconfigure the host_parameters array with the new indexing scheme.
    GOOFIT_TRACE("host_parameters[{}] = {}", host_parameters.size(), parametersList.size());
    parametersIdx = host_parameters.size();
    host_parameters.push_back(parametersList.size());
    for(auto &i : parametersList) {
        GOOFIT_TRACE("host_parameters[{}] = {}", host_parameters.size(), i.getValue());
        host_parameters.push_back(i.getValue());
    }

    GOOFIT_TRACE("host_constants[{}] = {}", host_constants.size(), constantsList.size());
    constantsIdx = host_constants.size();
    host_constants.push_back(constantsList.size());
    for(double i : constantsList) {
        GOOFIT_TRACE("host_constants[{}] = {}", host_constants.size(), i);
        host_constants.push_back(i);
    }

    GOOFIT_TRACE("host_observables[{}] = {}", host_observables.size(), observablesList.size());
    observablesIdx = host_observables.size();
    host_observables.push_back(observablesList.size());
    for(auto &i : observablesList) {
        GOOFIT_TRACE("host_observables[{}] = {}", host_observables.size(), i.getIndex());
        host_observables.push_back(i.getIndex());
    }

    GOOFIT_TRACE("host_normalizations[{}] = {}", host_normalizations.size(), 1);
    normalIdx = host_normalizations.size();
    host_normalizations.push_back(1.0);
    GOOFIT_TRACE("host_normalizations[{}] = {}", host_normalizations.size(), 0);
    host_normalizations.push_back(0.0);

    int numResonances = decayInfo.resonances.size();

    // add our resonance functions
    for(unsigned int i = 0; i < numResonances; i++)
        components[i]->recursiveSetIndices();

    // TODO: Add resolution function here
    resolutionFunction = host_function_table.size();
    components[numResonances]->recursiveSetIndices();

    // Next index starts our efficiency function
    efficiencyFunction = host_function_table.size();
    for(unsigned int i = numResonances + 1; i < components.size(); i++)
        components[i]->recursiveSetIndices();

    // update constants
    constantsList[constantsList.size() - 1] = host_function_table.size();
    GOOFIT_TRACE("Rewriting constants!");
    for(int i = 0; i < constantsList.size(); i++) {
        GOOFIT_TRACE("host_constants[{}] = {}", constantsIdx, constantsList[i]);
        host_constants[constantsIdx + 1 + i] = constantsList[i];
    }
    // TODO: This might be easy to clean up with the smart vectors.
}
__host__ void Amp3Body_TD::setDataSize(unsigned int dataSize, unsigned int evtSize, unsigned int offset) {
    // Default 5 is m12, m13, time, sigma_t, evtNum
    totalEventSize = evtSize;
    if(totalEventSize < 5)
        throw GooFit::GeneralError("totalEventSize {} must be 5 or more", totalEventSize);

    if(cachedWaves[0]) {
        for(auto &cachedWave : cachedWaves)
            delete cachedWave;
    }

    numEntries  = dataSize;
    eventOffset = offset;

// Ideally this would not be required, this would be called AFTER setData which will set m_iEventsPerTask
#ifdef GOOFIT_MPI
    int myId, numProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myId);

    int perTask = numEntries / numProcs;

    int *counts = new int[numProcs];

    for(int i = 0; i < numProcs - 1; i++)
        counts[i] = perTask;

    counts[numProcs - 1] = numEntries - perTask * (numProcs - 1);

    setNumPerTask(this, counts[myId]);

    delete[] counts;
#endif

    for(int i = 0; i < 16; i++) {
#ifdef GOOFIT_MPI
        cachedWaves[i] = new thrust::device_vector<WaveHolder_s>(m_iEventsPerTask);
#else
        cachedWaves[i] = new thrust::device_vector<WaveHolder_s>(dataSize);
#endif
        void *dummy = thrust::raw_pointer_cast(cachedWaves[i]->data());
        MEMCPY_TO_SYMBOL(cWaves,
                         &dummy,
                         sizeof(WaveHolder_s *),
                         ((16 * cacheToUse) + i) * sizeof(WaveHolder_s *),
                         cudaMemcpyHostToDevice);
    }

    setForceIntegrals();
}

__host__ void Amp3Body_TD::setD0Fraction(fptype d0fraction) {
    _D0Fraction = d0fraction;
    assert(_D0Fraction >= 0);
    assert(_D0Fraction <= 1);
}

__host__ fptype Amp3Body_TD::getD0Fraction() { return _D0Fraction; }

__host__ auto Amp3Body_TD::normalize() -> fptype {
    recursiveSetNormalization(1.0); // Not going to normalize efficiency,
    // so set normalization factor to 1 so it doesn't get multiplied by zero.
    // Copy at this time to ensure that the SpecialWaveCalculators, which need the efficiency,
    // don't get zeroes through multiplying by the normFactor.

    host_normalizations.sync(d_normalizations);

    int totalBins = _m12.getNumBins() * _m13.getNumBins();

    if(!dalitzNormRange) {
        gooMalloc((void **)&dalitzNormRange, 6 * sizeof(fptype));

        auto *host_norms = new fptype[6];
        host_norms[0]    = _m12.getLowerLimit();
        host_norms[1]    = _m12.getUpperLimit();
        host_norms[2]    = _m12.getNumBins();
        host_norms[3]    = _m13.getLowerLimit();
        host_norms[4]    = _m13.getUpperLimit();
        host_norms[5]    = _m13.getNumBins();
        MEMCPY(dalitzNormRange, host_norms, 6 * sizeof(fptype), cudaMemcpyHostToDevice);
        delete[] host_norms;
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

    // NB, SpecialWaveCalculator assumes that fit is unbinned!
    // And it needs to know the total event size, not just observables
    // for this particular PDF component.
    thrust::constant_iterator<fptype *> dataArray(dev_event_array);
    thrust::constant_iterator<int> eventSize(totalEventSize);
    thrust::counting_iterator<int> eventIndex(eventOffset);

    static int normCall = 0;
    normCall++;

    for(int i = 0; i < decayInfo.resonances.size(); ++i) {
        // printf("calculate i=%i, res_i=%i\n", i, decayInfo->resonances[i]->getFunctionIndex());
        calculators[i]->setTddpIndex(getFunctionIndex());
        calculators[i]->setResonanceIndex(decayInfo.resonances[i]->getFunctionIndex());
        if(redoIntegral[i]) {
#ifdef GOOFIT_MPI
            thrust::transform(
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex + m_iEventsPerTask, dataArray, eventSize)),
                strided_range<thrust::device_vector<WaveHolder_s>::iterator>(
                    cachedWaves[i]->begin(), cachedWaves[i]->end(), 1)
                    .begin(),
                *(calculators[i]));
#else
            thrust::transform(
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, dataArray, eventSize)),
                strided_range<thrust::device_vector<WaveHolder_s>::iterator>(
                    cachedWaves[i]->begin(), cachedWaves[i]->end(), 1)
                    .begin(),
                *(calculators[i]));
#endif
        }

        // Possibly this can be done more efficiently by exploiting symmetry?
        for(int j = 0; j < decayInfo.resonances.size(); ++j) {
            if((!redoIntegral[i]) && (!redoIntegral[j]))
                continue;

            integrators[i][j]->setTddpIndex(getFunctionIndex());
            integrators[i][j]->setResonanceIndex(decayInfo.resonances[i]->getFunctionIndex());
            integrators[i][j]->setEfficiencyIndex(decayInfo.resonances[j]->getFunctionIndex());

            ThreeComplex dummy(0, 0, 0, 0, 0, 0);
            SpecialComplexSum complexSum;
            thrust::constant_iterator<int> effFunc(efficiencyFunction);
            (*(integrals[i][j])) = thrust::transform_reduce(
                thrust::make_zip_iterator(thrust::make_tuple(binIndex, arrayAddress, effFunc)),
                thrust::make_zip_iterator(thrust::make_tuple(binIndex + totalBins, arrayAddress, effFunc)),
                *(integrators[i][j]),
                dummy,
                complexSum);
        }
    }

    // End of time-consuming integrals.

    fpcomplex integralA_2(0, 0);
    fpcomplex integralB_2(0, 0);
    fpcomplex integralABs(0, 0);

    for(unsigned int i = 0; i < decayInfo.resonances.size(); ++i) {
        fpcomplex amplitude_i(host_parameters[parametersIdx + i * 2 + 6], host_parameters[parametersIdx + i * 2 + 7]);

        for(unsigned int j = 0; j < decayInfo.resonances.size(); ++j) {
            fpcomplex amplitude_j(host_parameters[parametersIdx + j * 2 + 6],
                                  -host_parameters[parametersIdx + j * 2 + 7]); // Notice complex conjugation

            integralA_2 += (amplitude_i * amplitude_j
                            * fpcomplex(thrust::get<0>(*(integrals[i][j])), thrust::get<1>(*(integrals[i][j]))));
            integralABs += (amplitude_i * amplitude_j
                            * fpcomplex(thrust::get<2>(*(integrals[i][j])), thrust::get<3>(*(integrals[i][j]))));
            integralB_2 += (amplitude_i * amplitude_j
                            * fpcomplex(thrust::get<4>(*(integrals[i][j])), thrust::get<5>(*(integrals[i][j]))));
        }
    }

    double dalitzIntegralOne = integralA_2.real(); // Notice that this is already the abs2, so it's real by
                                                   // construction; but the compiler doesn't know that.
    double dalitzIntegralTwo = integralB_2.real();
    double dalitzIntegralThr = integralABs.real();
    double dalitzIntegralFou = integralABs.imag();

    fptype tau           = host_parameters[parametersIdx + 1];
    fptype xmixing0      = host_parameters[parametersIdx + 2];
    fptype ymixing0      = host_parameters[parametersIdx + 3];
    fptype deltax        = host_parameters[parametersIdx + 4];
    fptype deltay        = host_parameters[parametersIdx + 5];
    fptype xmixing_D0    = xmixing0 + deltax;
    fptype ymixing_D0    = ymixing0 + deltay;
    fptype xmixing_D0bar = xmixing0 - deltax;
    fptype ymixing_D0bar = ymixing0 - deltay;

    fptype ret_D0 = resolution->normalization(
        dalitzIntegralOne, dalitzIntegralTwo, dalitzIntegralThr, dalitzIntegralFou, tau, xmixing_D0, ymixing_D0);

    fptype ret_D0bar = resolution->normalization(
        dalitzIntegralOne, dalitzIntegralTwo, dalitzIntegralThr, dalitzIntegralFou, tau, xmixing_D0bar, ymixing_D0bar);

    // fptype _D0Fraction = 0.5; // Set D0 fraction to 1 for now.
    fptype ret = _D0Fraction * ret_D0 + (1. - _D0Fraction) * ret_D0bar;

    fptype xfix   = 0.0039;
    fptype yfix   = 0.0065;
    fptype taufix = 0.4101;

    if(_mistag) {
        ret *= (1 - 2 * _mistag.getValue());
        ret += _mistag.getValue()
               * resolution->normalization(
                   dalitzIntegralOne, dalitzIntegralTwo, dalitzIntegralThr, dalitzIntegralFou, taufix, xfix, yfix);
        ret += _mistag.getValue()
               * resolution->normalization(
                   dalitzIntegralTwo, dalitzIntegralOne, dalitzIntegralThr, -dalitzIntegralFou, taufix, xfix, yfix);
    }

    double binSizeFactor = 1;
    binSizeFactor *= ((_m12.getUpperLimit() - _m12.getLowerLimit()) / _m12.getNumBins());
    binSizeFactor *= ((_m13.getUpperLimit() - _m13.getLowerLimit()) / _m13.getNumBins());
    ret *= binSizeFactor;

    host_normalizations.at(normalIdx + 1) = 1.0 / ret;
    cachedNormalization                   = 1.0 / ret;

    GOOFIT_TRACE("Normalisation: {}", ret);

    return ret;
}

__host__ std::vector<std::vector<fptype>> Amp3Body_TD::getFractions() {
    recursiveSetNormalization(1.0); // Not going to normalize efficiency,
    // so set normalization factor to 1 so it doesn't get multiplied by zero.
    // Copy at this time to ensure that the SpecialWaveCalculators, which need the efficiency,
    // don't get zeroes through multiplying by the normFactor.

    host_normalizations.sync(d_normalizations);

    int totalBins = _m12.getNumBins() * _m13.getNumBins();

    if(!dalitzNormRange) {
        gooMalloc((void **)&dalitzNormRange, 6 * sizeof(fptype));

        auto *host_norms = new fptype[6];
        host_norms[0]    = _m12.getLowerLimit();
        host_norms[1]    = _m12.getUpperLimit();
        host_norms[2]    = _m12.getNumBins();
        host_norms[3]    = _m13.getLowerLimit();
        host_norms[4]    = _m13.getUpperLimit();
        host_norms[5]    = _m13.getNumBins();
        MEMCPY(dalitzNormRange, host_norms, 6 * sizeof(fptype), cudaMemcpyHostToDevice);
        delete[] host_norms;
    }

    std::vector<fptype> fracLists;
    fracLists.clear();

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

    // NB, SpecialWaveCalculator assumes that fit is unbinned!
    // And it needs to know the total event size, not just observables
    // for this particular PDF component.
    thrust::constant_iterator<fptype *> dataArray(dev_event_array);
    thrust::constant_iterator<int> eventSize(totalEventSize);
    thrust::counting_iterator<int> eventIndex(0);

    static int normCall = 0;
    normCall++;

    for(int i = 0; i < decayInfo.resonances.size(); ++i) {
        // printf("calculate i=%i, res_i=%i\n", i, decayInfo->resonances[i]->getFunctionIndex());
        calculators[i]->setTddpIndex(getFunctionIndex());
        calculators[i]->setResonanceIndex(decayInfo.resonances[i]->getFunctionIndex());
        if(redoIntegral[i]) {
#ifdef GOOFIT_MPI
            thrust::transform(
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex + m_iEventsPerTask, arrayAddress, eventSize)),
                strided_range<thrust::device_vector<WaveHolder_s>::iterator>(
                    cachedWaves[i]->begin(), cachedWaves[i]->end(), 1)
                    .begin(),
                *(calculators[i]));
#else
            thrust::transform(
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, arrayAddress, eventSize)),
                strided_range<thrust::device_vector<WaveHolder_s>::iterator>(
                    cachedWaves[i]->begin(), cachedWaves[i]->end(), 1)
                    .begin(),
                *(calculators[i]));
#endif
        }
        // Possibly this can be done more efficiently by exploiting symmetry?
        for(int j = 0; j < decayInfo.resonances.size(); ++j) {
            if((!redoIntegral[i]) && (!redoIntegral[j]))
                continue;

            integrators[i][j]->setTddpIndex(getFunctionIndex());
            integrators[i][j]->setResonanceIndex(decayInfo.resonances[i]->getFunctionIndex());
            integrators[i][j]->setEfficiencyIndex(decayInfo.resonances[j]->getFunctionIndex());

            ThreeComplex dummy(0, 0, 0, 0, 0, 0);
            SpecialComplexSum complexSum;
            thrust::constant_iterator<int> effFunc(efficiencyFunction);
            (*(integrals[i][j])) = thrust::transform_reduce(
                thrust::make_zip_iterator(thrust::make_tuple(binIndex, arrayAddress, effFunc)),
                thrust::make_zip_iterator(thrust::make_tuple(binIndex + totalBins, arrayAddress, effFunc)),
                *(integrators[i][j]),
                dummy,
                complexSum);
            if(i != j)
                (*(integrals[j][i])) = ThreeComplex(thrust::get<0>(*(integrals[i][j])),
                                                    -thrust::get<1>(*(integrals[i][j])),
                                                    thrust::get<2>(*(integrals[i][j])),
                                                    -thrust::get<3>(*(integrals[i][j])),
                                                    thrust::get<4>(*(integrals[i][j])),
                                                    -thrust::get<5>(*(integrals[i][j])));
        }
    }

    // End of time-consuming integrals.

    fpcomplex integralA_2(0, 0);
    const unsigned int nres = decayInfo.resonances.size();
    fptype matdiag[nres];
    //    fptype matdiag_int[nres][nres];
    std::vector<std::vector<fptype>> matdiag_int(nres, std::vector<fptype>(nres));

    for(unsigned int i = 0; i < decayInfo.resonances.size(); ++i) {
        fpcomplex amplitude_i(host_parameters[parametersIdx + i * 2 + 6], host_parameters[parametersIdx + i * 2 + 7]);
        std::string resname = decayInfo.resonances[i]->getName();

        for(unsigned int j = 0; j < decayInfo.resonances.size(); ++j) {
            fpcomplex amplitude_j(host_parameters[parametersIdx + j * 2 + 6],
                                  -host_parameters[parametersIdx + j * 2 + 7]); // Notice complex conjugation

            integralA_2 += (amplitude_i * amplitude_j
                            * fpcomplex(thrust::get<0>(*(integrals[i][j])), thrust::get<1>(*(integrals[i][j]))));
            if(i == j)
                matdiag[i] = (amplitude_i * amplitude_j
                              * fpcomplex(thrust::get<0>(*(integrals[i][j])), thrust::get<1>(*(integrals[i][j]))))
                                 .real();
            matdiag_int[i][j] = (amplitude_i * amplitude_j
                                 * fpcomplex(thrust::get<0>(*(integrals[i][j])),
                                             thrust::get<1>(*(integrals[i][j]))))
                                    .real(); // for reporting interference fractions
        }
    }

    std::streamsize ss = std::cout.precision();

    for(unsigned int i = 0; i < decayInfo.resonances.size(); ++i) {
        std::cout << std::setprecision(4) << "Integral contribution for res # " << i << " ( "
                  << decayInfo.resonances[i]->getName() << ") : " << (matdiag[i] / (integralA_2).real()) * 100. << "%."
                  << std::endl;
        fracLists.push_back(matdiag[i] / (integralA_2).real());
    }

    // MW 2 Aug 2016. Add interference fractions.
    for(unsigned int i = 0; i < decayInfo.resonances.size(); ++i) {
        for(unsigned int j = i + 1; j < decayInfo.resonances.size(); ++j) {
            std::cout << std::setprecision(4) << "Integral interference fraction for resonances " << i << ", " << j
                      << " ( " << decayInfo.resonances[i]->getName() << ", " << decayInfo.resonances[j]->getName()
                      << ") : " << (matdiag_int[i][j] / (integralA_2).real()) * 100. << "%." << std::endl;
        }
    }
    std::cout.precision(ss);

    // return (integralA_2).real();
    return matdiag_int;
}

} // namespace GooFit
