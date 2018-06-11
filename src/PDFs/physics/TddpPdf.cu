#include <goofit/Error.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/TddpPdf.h>

#include <thrust/transform_reduce.h>

#ifdef GOOFIT_MPI
#include <mpi.h>
#endif

namespace GooFit {

const int resonanceOffset = 8; // Offset of the first resonance into the parameter index array
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
__device__ WaveHolder_s *cWaves[16];

__device__ inline int parIndexFromResIndex(int resIndex) { return resonanceOffset + resIndex * resonanceSize; }

__device__ fpcomplex getResonanceAmplitude(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) {
    auto func = reinterpret_cast<resonance_function_ptr>(device_function_table[pc.funcIdx]);
    return (*func)(m12, m13, m23, pc);
}

__device__ ThreeComplex
device_Tddp_calcIntegrals(fptype m12, fptype m13, int res_i, int res_j, ParameterContainer &pc) {
    // For calculating Dalitz-plot integrals. What's needed is the products
    // AiAj*, AiBj*, and BiBj*, where
    // Ai = BW_i(x, y) + BW_i(y, x)
    // and Bi reverses the sign of the second BW.
    // This function returns the above values at a single point.
    // NB: Multiplication by efficiency is done by the calling function.
    // Note that this function expects
    // to be called on a normalization grid, not on
    // observed points, that's why it doesn't use
    // cWaves. No need to cache the values at individual
    // grid points - we only care about totals.

    ThreeComplex ret;

    if(!inDalitz(m12, m13, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass))
        return ret;

    fptype m23 = c_motherMass * c_motherMass + c_daug1Mass * c_daug1Mass + c_daug2Mass * c_daug2Mass
                 + c_daug3Mass * c_daug3Mass - m12 - m13;

    ParameterContainer ipc = pc;
    while(ipc.funcIdx < res_i)
        ipc.incrementIndex();

    ParameterContainer t = ipc;
    fpcomplex ai         = getResonanceAmplitude(m12, m13, m23, t);
    t                    = ipc;
    fpcomplex bi         = getResonanceAmplitude(m13, m12, m23, t);

    ParameterContainer jpc = pc;
    while(jpc.funcIdx < res_j)
        jpc.incrementIndex();

    t            = jpc;
    fpcomplex aj = conj(getResonanceAmplitude(m12, m13, m23, t));
    t            = jpc;
    fpcomplex bj = conj(getResonanceAmplitude(m13, m12, m23, t));

    ret = ThreeComplex(
        (ai * aj).real(), (ai * aj).imag(), (ai * bj).real(), (ai * bj).imag(), (bi * bj).real(), (bi * bj).imag());
    return ret;
}

__device__ fptype device_Tddp(fptype *evt, ParameterContainer &pc) {
    int num_parameters  = pc.getNumParameters();
    int num_constants   = pc.getNumConstants();
    int num_observables = pc.getNumObservables();

    int id_m12 = pc.getObservable(2);
    int id_m13 = pc.getObservable(3);
    int id_num = pc.getObservable(4);
    int id_mis = 0;
    if(num_observables > 5)
        id_mis = pc.getObservable(5);

    fptype m12 = evt[id_m12];
    fptype m13 = evt[id_m13];

    unsigned int numResonances = pc.getConstant(0);

    if(!inDalitz(m12, m13, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass)) {
        unsigned int endEfficiencyFunc = RO_CACHE(pc.constants[pc.constantIdx + 4]);
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

    auto evtNum = static_cast<int>(floor(0.5 + evt[id_num]));

    fpcomplex sumWavesA(0, 0);
    fpcomplex sumWavesB(0, 0);
    fpcomplex sumRateAA(0, 0);
    fpcomplex sumRateAB(0, 0);
    fpcomplex sumRateBB(0, 0);

    // unsigned int cacheToUse = pc.getConstant(1);
    fptype mistag = pc.getConstant(2);

    for(int i = 0; i < numResonances; ++i) {
        // int paramIndex = parIndexFromResIndex(i);
        fpcomplex amp{pc.getParameter(i * 2 + 3), pc.getParameter(i * 2 + 4)};

        // fpcomplex matrixelement(thrust::get<0>(cWaves[cacheToUse][evtNum*numResonances + i]),
        //				     thrust::get<1>(cWaves[cacheToUse][evtNum*numResonances + i]));
        // Note, to make this more efficient we should change it to only an array of fptype's, and read double2 at a
        // time.
        fpcomplex ai{RO_CACHE(cWaves[i][evtNum].ai_real), RO_CACHE(cWaves[i][evtNum].ai_imag)};
        fpcomplex bi{RO_CACHE(cWaves[i][evtNum].bi_real), RO_CACHE(cWaves[i][evtNum].bi_imag)};

        fpcomplex matrixelement = ai * amp;
        sumWavesA += matrixelement;

        // matrixelement = fpcomplex(thrust::get<2>(cWaves[cacheToUse][evtNum*numResonances + i]),
        //				       thrust::get<3>(cWaves[cacheToUse][evtNum*numResonances + i]));
        matrixelement = bi * amp;
        sumWavesB += matrixelement;
    }

    fptype _tau     = pc.getParameter(0);
    fptype _xmixing = pc.getParameter(1);
    fptype _ymixing = pc.getParameter(2);

    int id_time  = pc.getObservable(0);
    int id_sigma = pc.getObservable(1);

    fptype _time  = evt[id_time];
    fptype _sigma = evt[id_sigma];

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

    ret = (*(reinterpret_cast<device_resfunction_ptr>(device_function_table[pc.funcIdx])))(
        term1, term2, sumWavesA.real(), sumWavesA.imag(), _tau, _time, _xmixing, _ymixing, _sigma, pc);

    // For the reversed (mistagged) fraction, we make the
    // interchange A <-> B. So term1 stays the same,
    // term2 changes sign, and AB* becomes BA*.
    // Efficiency remains the same for the mistagged part,
    // because it depends on the momenta of the pi+ and pi-,
    // which don't change even though we tagged a D0 as D0bar.

    // fptype mistag = RO_CACHE(functorConstants[RO_CACHE(indices[1]) + 5]);

    if(mistag > 0) { // This should be either true or false for all events, so no branch is caused.
        // See header file for explanation of 'mistag' variable - it is actually the probability
        // of having the correct sign, given that we have a correctly reconstructed D meson.
        mistag = evt[id_mis];
        ret *= mistag;
        ret += (1 - mistag)
               * (*(reinterpret_cast<device_resfunction_ptr>(device_function_table[pc.funcIdx])))(
                     term1, -term2, sumWavesA.real(), -sumWavesA.imag(), _tau, _time, _xmixing, _ymixing, _sigma, pc);
    }

    // increment our resolution function
    pc.incrementIndex();

    fptype eff = callFunction(evt, pc);
    // internalDebug = 0;
    ret *= eff;

    return ret;
}

__device__ device_function_ptr ptr_to_Tddp = device_Tddp;

__host__ TddpPdf::TddpPdf(std::string n,
                          Observable _dtime,
                          Observable _sigmat,
                          Observable m12,
                          Observable m13,
                          EventNumber eventNumber,
                          DecayInfo3t decay,
                          MixingTimeResolution *r,
                          GooPdf *efficiency,
                          Observable *mistag)
    : GooPdf(n, _dtime, _sigmat, m12, m13, eventNumber)
    , decayInfo(decay)
    , _m12(m12)
    , _m13(m13)
    , resolution(r)
    , totalEventSize(6) // Default 5 = m12, m13, time, sigma_t, evtNum
{
    for(auto &cachedWave : cachedWaves)
        cachedWave = nullptr;

    if(mistag) {
        registerObservable(*mistag);
        totalEventSize = 6;
    }

    MEMCPY_TO_SYMBOL(c_motherMass, &decay.motherMass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug1Mass, &decay.daug1Mass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug2Mass, &decay.daug2Mass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug3Mass, &decay.daug3Mass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_meson_radius, &decay.meson_radius, sizeof(fptype), 0, cudaMemcpyHostToDevice);

    registerParameter(decay._tau);
    registerParameter(decay._xmixing);
    registerParameter(decay._ymixing);

    if(resolution->getDeviceFunction() < 0)
        throw GooFit::GeneralError("The resolution device function index {} must be more than 0",
                                   resolution->getDeviceFunction());

    registerConstant(decay.resonances.size());

    static int cacheCount = 0;
    cacheToUse            = cacheCount++;
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

    // TODO: Figure out what this needs?
    resolution->createParameters(this);

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

    addSpecialMask(PdfBase::ForceSeparateNorm);
}

__host__ TddpPdf::TddpPdf(std::string n,
                          Observable _dtime,
                          Observable _sigmat,
                          Observable m12,
                          Observable m13,
                          EventNumber eventNumber,
                          DecayInfo3t decay,
                          std::vector<MixingTimeResolution *> &r,
                          GooPdf *efficiency,
                          Observable md0,
                          Observable *mistag)
    : GooPdf(n, _dtime, _sigmat, m12, m13, eventNumber, md0)
    , decayInfo(decay)
    , _m12(m12)
    , _m13(m13)
    , resolution(
          r[0]) // Only used for normalization, which only depends on x and y - it doesn't matter which one we use.
    , totalEventSize(6) // This case adds the D0 mass by default.
{
    for(auto &cachedWave : cachedWaves)
        cachedWave = nullptr;

    if(mistag) {
        registerObservable(*mistag);
        totalEventSize++;
    }

    MEMCPY_TO_SYMBOL(c_motherMass, &decay.motherMass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug1Mass, &decay.daug1Mass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug2Mass, &decay.daug2Mass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug3Mass, &decay.daug3Mass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_meson_radius, &decay.meson_radius, sizeof(fptype), 0, cudaMemcpyHostToDevice);

    registerParameter(decay._tau);
    registerParameter(decay._xmixing);
    registerParameter(decay._ymixing);
    printf("Multiple resolution functions not supported yet!\n");

    registerConstant(decayInfo.resonances.size());

    static int cacheCount = 0;
    cacheToUse            = cacheCount++;
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

        // TODO:
        i->createParameters(this);

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

    addSpecialMask(PdfBase::ForceSeparateNorm);
}

// Note: We need to manually populate the arrays so we can track the efficiency function!
__host__ void TddpPdf::populateArrays() {
    // populate all the arrays
    GOOFIT_TRACE("TddpPdf: Populating Arrays for {}", getName());

    // reconfigure the host_parameters array with the new indexing scheme.
    GOOFIT_TRACE("host_parameters[{}] = {}", totalParameters, parametersList.size());
    host_parameters[totalParameters] = parametersList.size();
    parametersIdx                    = totalParameters;
    totalParameters++;
    for(auto &i : parametersList) {
        GOOFIT_TRACE("host_parameters[{}] = {}", totalParameters, i.getValue());
        host_parameters[totalParameters] = i.getValue();
        totalParameters++;
    }

    GOOFIT_TRACE("host_constants[{}] = {}", totalConstants, constantsList.size());
    host_constants[totalConstants] = constantsList.size();
    constantsIdx                   = totalConstants;
    totalConstants++;
    for(double i : constantsList) {
        GOOFIT_TRACE("host_constants[{}] = {}", totalConstants, i);
        host_constants[totalConstants] = i;
        totalConstants++;
    }

    GOOFIT_TRACE("host_observables[{}] = {}", totalObservables, observablesList.size());
    host_observables[totalObservables] = observablesList.size();
    observablesIdx                     = totalObservables;
    totalObservables++;
    for(auto &i : observablesList) {
        GOOFIT_TRACE("host_observables[{}] = {}", totalObservables, i.getIndex());
        host_observables[totalObservables] = i.getIndex();
        totalObservables++;
    }

    GOOFIT_TRACE("host_normalizations[{}] = {}", totalNormalizations, 1);
    host_normalizations[totalNormalizations] = 1;
    normalIdx                                = totalNormalizations++;
    GOOFIT_TRACE("host_normalizations[{}] = {}", totalNormalizations, 0);
    host_normalizations[totalNormalizations] = cachedNormalization;
    totalNormalizations++;

    int numResonances = decayInfo.resonances.size();

    // add our resonance functions
    for(unsigned int i = 0; i < numResonances; i++)
        components[i]->recursiveSetIndices();

    // TODO: Add resolution function here
    resolutionFunction = num_device_functions;
    components[numResonances]->recursiveSetIndices();

    // Next index starts our efficiency function
    efficiencyFunction = num_device_functions;
    for(unsigned int i = numResonances + 1; i < components.size(); i++)
        components[i]->recursiveSetIndices();

    // update constants
    constantsList[constantsList.size() - 1] = num_device_functions;
    GOOFIT_TRACE("Rewriting constants!");
    for(int i = 0; i < constantsList.size(); i++) {
        GOOFIT_TRACE("host_constants[{}] = {}", constantsIdx, constantsList[i]);
        host_constants[constantsIdx + 1 + i] = constantsList[i];
    }
}
__host__ void TddpPdf::setDataSize(unsigned int dataSize, unsigned int evtSize) {
    // Default 5 is m12, m13, time, sigma_t, evtNum
    totalEventSize = evtSize;
    if(totalEventSize < 5)
        throw GooFit::GeneralError("totalEventSize {} must be 5 or more", totalEventSize);

    if(cachedWaves[0]) {
        for(auto &cachedWave : cachedWaves)
            delete cachedWave;
    }

    numEntries = dataSize;

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
        MEMCPY_TO_SYMBOL(cWaves, &dummy, sizeof(WaveHolder_s *), i * sizeof(WaveHolder_s *), cudaMemcpyHostToDevice);
    }

    setForceIntegrals();
}

__host__ fptype TddpPdf::normalize() {
    recursiveSetNormalization(1.0); // Not going to normalize efficiency,
    // so set normalization factor to 1 so it doesn't get multiplied by zero.
    // Copy at this time to ensure that the SpecialWaveCalculators, which need the efficiency,
    // don't get zeroes through multiplying by the normFactor.
    MEMCPY_TO_SYMBOL(
        d_normalizations, host_normalizations, totalNormalizations * sizeof(fptype), 0, cudaMemcpyHostToDevice);

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
            // std::cout << "Integral for resonance " << i << " " << numEntries << " " << totalEventSize << std::endl;
        }

        // Possibly this can be done more efficiently by exploiting symmetry?
        for(int j = 0; j < decayInfo.resonances.size(); ++j) {
            if((!redoIntegral[i]) && (!redoIntegral[j]))
                continue;

            integrators[i][j]->setTddpIndex(getFunctionIndex());
            integrators[i][j]->setResonanceIndex(decayInfo.resonances[i]->getFunctionIndex());
            integrators[i][j]->setEfficiencyIndex(decayInfo.resonances[j]->getFunctionIndex());

            // printf("integrate i=%i j=%i, res_i=%i res_j=%i\n", i, j, decayInfo->resonances[i]->getFunctionIndex (),
            //    decayInfo->resonances[j]->getFunctionIndex ());
            ThreeComplex dummy(0, 0, 0, 0, 0, 0);
            SpecialComplexSum complexSum;
            thrust::constant_iterator<int> effFunc(efficiencyFunction);
            (*(integrals[i][j])) = thrust::transform_reduce(
                thrust::make_zip_iterator(thrust::make_tuple(binIndex, arrayAddress, effFunc)),
                thrust::make_zip_iterator(thrust::make_tuple(binIndex + totalBins, arrayAddress, effFunc)),
                *(integrators[i][j]),
                dummy,
                complexSum);
            /*
            std::cout << "With resonance " << j << ": "
            << thrust::get<0>(*(integrals[i][j])) << " "
            << thrust::get<1>(*(integrals[i][j])) << " "
            << thrust::get<2>(*(integrals[i][j])) << " "
            << thrust::get<3>(*(integrals[i][j])) << " "
            << thrust::get<4>(*(integrals[i][j])) << " "
            << thrust::get<5>(*(integrals[i][j])) << std::endl;
            */
        }
    }

    // End of time-consuming integrals.

    fpcomplex integralA_2(0, 0);
    fpcomplex integralB_2(0, 0);
    fpcomplex integralABs(0, 0);

    for(unsigned int i = 0; i < decayInfo.resonances.size(); ++i) {
        // int param_i = parameters + resonanceOffset + resonanceSize * i;
        fpcomplex amplitude_i(host_parameters[parametersIdx + i * 2 + 4], host_parameters[parametersIdx + i * 2 + 5]);

        for(unsigned int j = 0; j < decayInfo.resonances.size(); ++j) {
            // int param_j = parameters + resonanceOffset + resonanceSize * j;
            fpcomplex amplitude_j(host_parameters[parametersIdx + j * 2 + 4],
                                  -host_parameters[parametersIdx + j * 2 + 5]); // Notice complex conjugation

            integralA_2 += (amplitude_i * amplitude_j
                            * fpcomplex(thrust::get<0>(*(integrals[i][j])), thrust::get<1>(*(integrals[i][j]))));
            integralABs += (amplitude_i * amplitude_j
                            * fpcomplex(thrust::get<2>(*(integrals[i][j])), thrust::get<3>(*(integrals[i][j]))));
            integralB_2 += (amplitude_i * amplitude_j
                            * fpcomplex(thrust::get<4>(*(integrals[i][j])), thrust::get<5>(*(integrals[i][j]))));

            /*
            if (cpuDebug & 1) {
            int idx = i * decayInfo.resonances.size() + j;
            if (0 == host_callnumber) std::cout << "Integral contribution " << i << ", " << j << " " << idx << " : "
                                << amplitude_i << " "
                                << amplitude_j << " ("
                                << real(amplitude_i * amplitude_j * complex<fptype>(thrust::get<0>(*(integrals[i][j])),
            thrust::get<1>(*(integrals[i][j])))) << ", "
                                << imag(amplitude_i * amplitude_j * complex<fptype>(thrust::get<0>(*(integrals[i][j])),
            thrust::get<1>(*(integrals[i][j])))) << ") ("
                                << real(amplitude_i * amplitude_j * complex<fptype>(thrust::get<2>(*(integrals[i][j])),
            thrust::get<3>(*(integrals[i][j])))) << ", "
                                << imag(amplitude_i * amplitude_j * complex<fptype>(thrust::get<2>(*(integrals[i][j])),
            thrust::get<3>(*(integrals[i][j])))) << ") ("
                                << real(amplitude_i * amplitude_j * complex<fptype>(thrust::get<4>(*(integrals[i][j])),
            thrust::get<5>(*(integrals[i][j])))) << ", "
                                << imag(amplitude_i * amplitude_j * complex<fptype>(thrust::get<4>(*(integrals[i][j])),
            thrust::get<5>(*(integrals[i][j])))) << ") "
                                << thrust::get<0>(*(integrals[i][j])) << ", "
                                << thrust::get<1>(*(integrals[i][j])) << ") ("
                                << thrust::get<2>(*(integrals[i][j])) << ", "
                                << thrust::get<3>(*(integrals[i][j])) << ") ("
                                << thrust::get<4>(*(integrals[i][j])) << ", "
                                << thrust::get<5>(*(integrals[i][j])) << ") ("
                                << real(integralA_2) << ", " << imag(integralA_2) << ") "
                                << std::endl;
                 }
                 */
        }
    }

    double dalitzIntegralOne = integralA_2.real(); // Notice that this is already the abs2, so it's real by
                                                   // construction; but the compiler doesn't know that.
    double dalitzIntegralTwo = integralB_2.real();
    double dalitzIntegralThr = integralABs.real();
    double dalitzIntegralFou = integralABs.imag();

    fptype tau     = host_parameters[parametersIdx + 1];
    fptype xmixing = host_parameters[parametersIdx + 2];
    fptype ymixing = host_parameters[parametersIdx + 3];

    fptype ret = resolution->normalization(
        dalitzIntegralOne, dalitzIntegralTwo, dalitzIntegralThr, dalitzIntegralFou, tau, xmixing, ymixing);

    double binSizeFactor = 1;
    binSizeFactor *= ((_m12.getUpperLimit() - _m12.getLowerLimit()) / _m12.getNumBins());
    binSizeFactor *= ((_m13.getUpperLimit() - _m13.getLowerLimit()) / _m13.getNumBins());
    ret *= binSizeFactor;

    host_normalizations[normalIdx + 1] = 1.0 / ret;
    cachedNormalization                = 1.0 / ret;
    // std::cout << "End of TDDP normalization: " << ret << " " << host_normalization[parameters] << " " <<
    // binSizeFactor << std::endl;
    return ret;
}
//#endif

SpecialDalitzIntegrator::SpecialDalitzIntegrator(int pIdx, unsigned int ri, unsigned int rj)
    : resonance_i(ri)
    , resonance_j(rj)
    , parameters(pIdx) {}

__device__ ThreeComplex SpecialDalitzIntegrator::operator()(thrust::tuple<int, fptype *, int> t) const {
    // Bin index, base address [lower, upper,getNumBins]
    // Notice that this is basically MetricTaker::operator (binned) with the special-case knowledge
    // that event size is two, and that the function to call is dev_Tddp_calcIntegrals.

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

    ParameterContainer pc;

    fptype events[10];

    // increment until we are at tddp index
    while(pc.funcIdx < tddp)
        pc.incrementIndex();

    int id_m12 = pc.getObservable(2);
    int id_m13 = pc.getObservable(3);
    // if (0 == THREADIDX) cuPrintf("%i %i %i %f %f operator\n", thrust::get<0>(t), thrust::get<0>(t) % numBinsM12,
    // globalBinNumber, binCenterM12, binCenterM13);
    ThreeComplex ret = device_Tddp_calcIntegrals(binCenterM12, binCenterM13, resonance_i, resonance_j, pc);

    // fptype fakeEvt[10]; // Need room for many observables in case m12 or m13 were assigned a high index in an
    // event-weighted fit.
    events[0]      = 2;
    events[id_m12] = binCenterM12;
    events[id_m13] = binCenterM13;
    // unsigned int numResonances                               = indices[6];
    // int effFunctionIdx                                       = parIndexFromResIndex(numResonances);
    // if (thrust::get<0>(t) == 19840) {internalDebug1 = BLOCKIDX; internalDebug2 = THREADIDX;}
    // fptype eff = (*(reinterpret_cast<device_function_ptr>(device_function_table[indices[effFunctionIdx]])))(fakeEvt,
    // cudaArray, paramIndices + indices[effFunctionIdx + 1]);
    while(pc.funcIdx < thrust::get<2>(t))
        pc.incrementIndex();

    fptype eff = callFunction(events, pc);
    // if (thrust::get<0>(t) == 19840) {
    // internalDebug1 = -1;
    // internalDebug2 = -1;
    // printf("Efficiency: %i %f %f %f %i\n", thrust::get<0>(t), binCenterM12, binCenterM13, eff, effFunctionIdx);
    // printf("Efficiency: %f %f %f %f %f %i %i\n", fakeEvt[0], fakeEvt[1], fakeEvt[2], fakeEvt[3], fakeEvt[4],
    // indices[indices[0] + 2 + 2], indices[indices[0] + 2 + 3]);
    //}

    // Multiplication by eff, not sqrt(eff), is correct:
    // These complex numbers will not be squared when they
    // go into the integrals. They've been squared already,
    // as it were.
    thrust::get<0>(ret) *= eff;
    thrust::get<1>(ret) *= eff;
    thrust::get<2>(ret) *= eff;
    thrust::get<3>(ret) *= eff;
    thrust::get<4>(ret) *= eff;
    thrust::get<5>(ret) *= eff;
    return ret;
}

SpecialWaveCalculator::SpecialWaveCalculator(int pIdx, unsigned int res_idx)
    : resonance_i(res_idx)
    , parameters(pIdx) {}

__device__ WaveHolder_s SpecialWaveCalculator::operator()(thrust::tuple<int, fptype *, int> t) const {
    // Calculates the BW values for a specific resonance.
    // The 'A' wave stores the value at each point, the 'B'
    // at the opposite (reversed) point.

    WaveHolder_s ret;
    ret.ai_real = 0.0;
    ret.ai_imag = 0.0;
    ret.bi_real = 0.0;
    ret.bi_imag = 0.0;

    int evtNum  = thrust::get<0>(t);
    int evtSize = thrust::get<2>(t);
    fptype *evt = thrust::get<1>(t) + (evtNum * evtSize);

    ParameterContainer pc;

    fptype events[10];

    for(int i = 0; i < evtSize; i++)
        events[i] = evt[i];

    // increment until we are at tddp index
    while(pc.funcIdx < tddp)
        pc.incrementIndex();

    int id_m12 = pc.getObservable(2);
    int id_m13 = pc.getObservable(3);

    // Read these values as tddp.
    fptype m12 = events[id_m12];
    fptype m13 = events[id_m13];

    if(!inDalitz(m12, m13, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass))
        return ret;

    fptype m23 = c_motherMass * c_motherMass + c_daug1Mass * c_daug1Mass + c_daug2Mass * c_daug2Mass
                 + c_daug3Mass * c_daug3Mass - m12 - m13;

    // int parameter_i       = parIndexFromResIndex(resonance_i); // Find position of this resonance relative to TDDP
    // start  unsigned int functn_i = indices[parameter_i + 2];  unsigned int params_i = indices[parameter_i + 3];

    while(pc.funcIdx < resonance_i)
        pc.incrementIndex();

    ParameterContainer tmp = pc;
    fpcomplex ai           = getResonanceAmplitude(m12, m13, m23, tmp);
    tmp                    = pc;
    fpcomplex bi           = getResonanceAmplitude(m13, m12, m23, tmp);

    // printf("Amplitudes %f, %f => (%f %f) (%f %f)\n", m12, m13, ai.real, ai.imag, bi.real, bi.imag);

    ret.ai_real = ai.real();
    ret.ai_imag = ai.imag();
    ret.bi_real = bi.real();
    ret.bi_imag = bi.imag();

    return ret;
}

} // namespace GooFit
