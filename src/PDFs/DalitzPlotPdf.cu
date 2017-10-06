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

__device__ thrust::complex<fptype> device_DalitzPlot_calcIntegrals(fptype m12, fptype m13, int res_i, int res_j, ParameterContainer &pc) {
    // Calculates BW_i(m12, m13) * BW_j^*(m12, m13).
    // This calculation is in a separate function so
    // it can be cached. Note that this function expects
    // to be called on a normalisation grid, not on
    // observed points, that's why it doesn't use
    // cResonances. No need to cache the values at individual
    // grid points - we only care about totals.
    fptype motherMass = c_motherMass;//RO_CACHE(pc.constants[pc.constantIdx + 4]);
    fptype daug1Mass  = c_daug1Mass;//RO_CACHE(pc.constants[pc.constantIdx + 5]);
    fptype daug2Mass  = c_daug2Mass;//RO_CACHE(pc.constants[pc.constantIdx + 6]);
    fptype daug3Mass  = c_daug3Mass;//RO_CACHE(pc.constants[pc.constantIdx + 7]);

    thrust::complex<fptype> ret;

    if(!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass))
        return ret;

    fptype m23 = motherMass * motherMass + daug1Mass * daug1Mass + daug2Mass * daug2Mass + daug3Mass * daug3Mass - m12 - m13;

    ParameterContainer ipc = pc;
    for (int i = 0; i < res_i; i++)
        ipc.incrementIndex ();
    //int parameter_i       = parIndexFromResIndex_DP(res_i);
    //unsigned int functn_i = RO_CACHE(indices[parameter_i + 2]);
    //unsigned int params_i = RO_CACHE(indices[parameter_i + 3]);
    
    ret                   = getResonanceAmplitude(m12, m13, m23, ipc);

    ParameterContainer jpc = pc;
    for (int i = 0; i < res_j; i++)
        jpc.incrementIndex ();
    //int parameter_j       = parIndexFromResIndex_DP(res_j);
    //unsigned int functn_j = RO_CACHE(indices[parameter_j + 2]);
    //unsigned int params_j = RO_CACHE(indices[parameter_j + 3]);
    ret *= conj(getResonanceAmplitude(m12, m13, m23, jpc));

    return ret;
}

__device__ fptype device_DalitzPlot(fptype *evt, ParameterContainer &pc) {
    int num_obs = RO_CACHE(pc.observables[pc.observableIdx]);
    int id_m12 = RO_CACHE(pc.observables[pc.observableIdx + 1]);
    int id_m13 = RO_CACHE(pc.observables[pc.observableIdx + 2]);
    int id_num = RO_CACHE(pc.observables[pc.observableIdx + 3]);

    //fptype motherMass = c_motherMass;//RO_CACHE(pc.constants[pc.constantIdx + 4]);
    //fptype daug1Mass  = c_daug1Mass;//RO_CACHE(pc.constants[pc.constantIdx + 5]);
    //fptype daug2Mass  = c_daug2Mass;//RO_CACHE(pc.constants[pc.constantIdx + 6]);
    //fptype daug3Mass  = c_daug3Mass;//RO_CACHE(pc.constants[pc.constantIdx + 7]);

    fptype m12 = RO_CACHE(evt[id_m12]);
    fptype m13 = RO_CACHE(evt[id_m13]);

    unsigned int numResonances = RO_CACHE(pc.constants[pc.constantIdx + 1]);
    unsigned int cacheToUse    = RO_CACHE(pc.constants[pc.constantIdx + 2]);

    if(!inDalitz(m12, m13, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass)) {
        pc.incrementIndex(1, numResonances*2, 2, num_obs, 1);

        //TODO: loop over resonances and efficiency functions
        return 0;
    }

    fptype evtIndex = RO_CACHE(evt[id_num]);

    auto evtNum = static_cast<int>(floor(0.5 + evtIndex));

    thrust::complex<fptype> totalAmp(0, 0);

    for(int i = 0; i < numResonances; ++i) {
        //int paramIndex              = parIndexFromResIndex_DP(i);
        //double2 *resPtr = reinterpret_cast<double2*> (cResonances[i]);
        thrust::complex<fptype> amp = thrust::complex<fptype>(RO_CACHE(pc.parameters[pc.parameterIdx + i*2 + 1]),
                                                              RO_CACHE(pc.parameters[pc.parameterIdx + i*2 + 2]));

        // potential performance improvement by
        //double2 *t = RO_CACHE(reinterpret_cast<double2*> (&(cResonances[i][evtNum])));
        //thrust::complex<fptype> me(t->x, t->y);
        //thrust::complex<fptype> me = RO_CACHE(cResonances[i][evtNum]);

        //double2 *ptr = reinterpret_cast<double2*> (cResonances[i][evtNum]);

        //double2 v = RO_CACHE(resPtr[evtNum]);

        //fptype me_real = cResonances[i][evtNum].real();
        //fptype me_imag = cResonances[i][evtNum].imag();
        //thrust::complex<fptype> me = cResonances[i][evtNum];
        //thrust::complex<fptype> me (me_real, me_imag);
        thrust::complex<fptype> me = RO_CACHE(cResonances[i][evtNum]);
        //thrust::complex<fptype> me (v.x, v.y);

        totalAmp += amp * me;
    }

    fptype ret         = thrust::norm(totalAmp);

    pc.incrementIndex(1, numResonances*2, 2, num_obs, 1);

    //loop to efficiency idx
    for (int i = 0; i < numResonances; i++)
        pc.incrementIndex ();

    fptype eff         = callFunction(evt, pc);
    ret *= eff;

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
    //MEMCPY_TO_SYMBOL(functorConstants, decayConstants, 5 * sizeof(fptype), cIndex * sizeof(fptype), cudaMemcpyHostToDevice);

    //Passing values to the defined constants.  Rather than push into list, which means each resonance
    //1. duplicates 5 fptypes
    //2. also has to index into the array
    //we instead just have a define.
    MEMCPY_TO_SYMBOL(c_motherMass, &decayConstants[0], sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug1Mass, &decayConstants[1], sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug2Mass, &decayConstants[2], sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug3Mass, &decayConstants[3], sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_meson_radius, &decayConstants[4], sizeof(fptype), 0, cudaMemcpyHostToDevice);

    pindices.push_back(decayInfo->resonances.size());
    constantsList.push_back (decayInfo->resonances.size());
    static int cacheCount = 0;
    cacheToUse            = cacheCount++;
    pindices.push_back(cacheToUse);
    constantsList.push_back(cacheToUse);

    for(auto &resonance : decayInfo->resonances) {
        pindices.push_back(registerParameter(resonance->amp_real));
        pindices.push_back(registerParameter(resonance->amp_imag));
        //pindices.push_back(resonance->getFunctionIndex());
        //pindices.push_back(resonance->getParameterIndex());
        //resonance->setConstantIndex(cIndex);
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

void DalitzPlotPdf::recursiveSetIndices () {
    GET_FUNCTION_ADDR(ptr_to_DalitzPlot);

    GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName (), "ptr_to_DalitzPlot");
    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx = num_device_functions++;
 
    populateArrays();

    //save our efficiency function.  Resonance's are saved first, then the efficiency function.  Take -1 as efficiency!
    efficiencyFunction = num_device_functions - 1;
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
    //we need to update the normal here, as values are used at this point.
    MEMCPY_TO_SYMBOL(d_normalisations, host_normalisations, totalNormalisations*sizeof(fptype), 0, cudaMemcpyHostToDevice);


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
	//grab the index for this resonance.
        calculators[i]->setResonanceIndex (decayInfo->resonances[i]->getFunctionIndex ());
        //calculators[i]->setEfficiencyIndex (efficiencyFunction);
        if(redoIntegral[i]) {
#ifdef GOOFIT_MPI
            thrust::transform(
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex + m_iEventsPerTask, arrayAddress, eventSize)),
                strided_range<thrust::device_vector<thrust::complex<fptype>>::iterator>(
                    cachedWaves[i]->begin(), cachedWaves[i]->end(), 1).begin(),
                *(calculators[i]));
#else
            thrust::transform(
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, arrayAddress, eventSize)),
                strided_range<thrust::device_vector<thrust::complex<fptype>>::iterator>(
                    cachedWaves[i]->begin(), cachedWaves[i]->end(), 1).begin(),
                *(calculators[i]));
#endif
        }

        // Possibly this can be done more efficiently by exploiting symmetry?
        for(int j = 0; j < decayInfo->resonances.size(); ++j) {
            if((!redoIntegral[i]) && (!redoIntegral[j]))
                continue;

            integrators[i][j]->setResonanceIndex(decayInfo->resonances[i]->getFunctionIndex());
            //integrators[i][j]->setEfficiencyIndex(efficiencyFunction);
            integrators[i][j]->setEfficiencyIndex(decayInfo->resonances[j]->getFunctionIndex());
            thrust::constant_iterator<int> effFunc(efficiencyFunction);

            thrust::complex<fptype> dummy(0, 0);
            thrust::plus<thrust::complex<fptype>> complexSum;
            (*(integrals[i][j])) = thrust::transform_reduce(
                thrust::make_zip_iterator(thrust::make_tuple(binIndex, arrayAddress, effFunc)),
                thrust::make_zip_iterator(thrust::make_tuple(binIndex + totalBins, arrayAddress, effFunc)),
                *(integrators[i][j]),
                dummy,
                complexSum);
        }
    }

    // End of time-consuming integrals.
    thrust::complex<fptype> sumIntegral(0, 0);

    for(unsigned int i = 0; i < decayInfo->resonances.size(); ++i) {
        //int param_i = parameters + resonanceOffset_DP + resonanceSize * i;
        thrust::complex<fptype> amplitude_i(host_parameters[parametersIdx + i*2 + 1], host_parameters[parametersIdx + i*2 + 2]);

	//printf("i:%i - %f,%f\n", i, amplitude_i.real(), amplitude_i.imag());

        for(unsigned int j = 0; j < decayInfo->resonances.size(); ++j) {
            //int param_j = parameters + resonanceOffset_DP + resonanceSize * j;
            thrust::complex<fptype> amplitude_j(host_parameters[parametersIdx + j*2 + 1], -host_parameters[parametersIdx + j*2 + 2]);

	    //printf("j:%i - %f,%f\n", j, amplitude_j.real(), amplitude_j.imag());
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

    host_normalisations[normalIdx + 1] = 1.0 / ret;
    // printf("%f %f\n", ret, binSizeFactor);
    return ret;
}

SpecialResonanceIntegrator::SpecialResonanceIntegrator(int pIdx, unsigned int ri, unsigned int rj)
    : resonance_i(ri)
    , resonance_j(rj)
    , parameters(pIdx) {}

__device__ thrust::complex<fptype> SpecialResonanceIntegrator::operator()(thrust::tuple<int, fptype *, int> t) const {
    //(brad): new indexing plan: bin number, function id, parameter id (not required), fptype with actual bins(needed???)

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
    auto numBinsM13      = static_cast<int>(floor(thrust::get<1>(t)[2] + 0.5));
    fptype binCenterM13  = upperBoundM13 - lowerBoundM13;
    binCenterM13 /= numBinsM13;
    binCenterM13 *= (globalBinNumber + 0.5);
    binCenterM13 += lowerBoundM13;

    ParameterContainer pc;

    thrust::complex<fptype> ret = device_DalitzPlot_calcIntegrals(binCenterM12, binCenterM13, resonance_i, resonance_j, pc);

    //TODO: read id's in in order to set them for the fake event.

    fptype fakeEvt[10]; // Need room for many observables in case m12 or m13 were assigned a high index in an
                        // event-weighted fit.
    fakeEvt[0] = 2;
    fakeEvt[1] = binCenterM12;
    fakeEvt[2] = binCenterM13;

    // unsigned int numResonances           = indices[2];
    // int effFunctionIdx                   = parIndexFromResIndex_DP(numResonances);

    //increment until we are on the efficiency function (17)
    while (pc.funcIdx < thrust::get<2>(t))
        pc.incrementIndex ();

    fptype eff = callFunction(fakeEvt, pc);

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

    fptype m12            = evt[0];
    fptype m13            = evt[1];

    //TODO: This function will need to find the base dalitz function
    ParameterContainer pc;

    fptype motherMass = c_motherMass;//pc.constants[pc.constantIdx + 4];
    fptype daug1Mass  = c_daug1Mass;//pc.constants[pc.constantIdx + 5];
    fptype daug2Mass  = c_daug2Mass;//pc.constants[pc.constantIdx + 6];
    fptype daug3Mass  = c_daug3Mass;//pc.constants[pc.constantIdx + 7];

    if(!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass))
        return ret;

    fptype m23 = motherMass * motherMass + daug1Mass * daug1Mass + daug2Mass * daug2Mass + daug3Mass * daug3Mass - m12 - m13;

    for (int i = 0; i < resonance_i; i++)
        pc.incrementIndex();

    ret = getResonanceAmplitude(m12, m13, m23, pc);

    return ret;
}

} // namespace GooFit
