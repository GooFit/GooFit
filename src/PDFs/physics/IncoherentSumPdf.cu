#include <goofit/Error.h>
#include <goofit/PDFs/physics/IncoherentSumPdf.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <thrust/complex.h>

#include <thrust/transform_reduce.h>

namespace GooFit {

const int resonanceOffset_incoherent = 4; // Offset of the first resonance into the parameter index array.
// Notice that this is different from the TddpPdf case because there's no time information.
// In particular the offset consists of nP, constant index, number of resonances, and cache index.

__device__ fpcomplex *cResonanceValues[10];

__device__ inline int parIndexFromResIndex_incoherent(int resIndex) {
    return resonanceOffset_incoherent + resIndex * resonanceSize;
}

__device__ fptype device_incoherent(fptype *evt, ParameterContainer &pc) {
    // Calculates the incoherent sum over the resonances.
    int numObs  = pc.getNumObservables();
    int evtId   = pc.getObservable(2);
    auto evtNum = static_cast<int>(floor(0.5 + evt[evtId]));

    fptype ret                 = 0;
    unsigned int numResonances = pc.getConstant(4);
    unsigned int cacheToUse    = pc.getConstant(5);

    for(int i = 0; i < numResonances; ++i) {
        // int paramIndex   = parIndexFromResIndex_incoherent(i);
        // fptype amplitude = p[indices[paramIndex + 0]];
        fptype amplitude = pc.getParameter(i);

        fpcomplex matrixelement = cResonanceValues[cacheToUse][evtNum * numResonances + i];
        ret += amplitude * thrust::norm(matrixelement);
    }

    // pc.incrementIndex(1, numResonances, 2, numObs, 1);
    pc.incrementIndex();

    // increment through resonances
    for(int i = 0; i < numResonances; i++)
        pc.incrementIndex();
    // Multiply by efficiency
    // int effFunctionIdx = parIndexFromResIndex_incoherent(numResonances);
    fptype eff = callFunction(evt, pc);

    ret *= eff;

    return ret;
}

__device__ device_function_ptr ptr_to_incoherent = device_incoherent;

__host__ IncoherentSumPdf::IncoherentSumPdf(
    std::string n, Observable m12, Observable m13, EventNumber eventNumber, DecayInfo3 decay, GooPdf *eff)
    : GooPdf(n, m12, m13, eventNumber)
    , decayInfo(decay)
    , _m12(m12)
    , _m13(m13)
    , dalitzNormRange(nullptr)
    , cachedResonances(nullptr)
    , integrals(nullptr)
    , forceRedoIntegrals(true)
    , totalEventSize(3) // Default 3 = m12, m13, evtNum. Will likely be overridden.
    , cacheToUse(0)
    , efficiency(eff)
    , integrators(nullptr)
    , calculators(nullptr) {
    registerConstant(observablesList.size());
    registerConstant(0);
    registerConstant(0);
    registerConstant(0);

    MEMCPY_TO_SYMBOL(c_motherMass, &decayInfo.motherMass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug1Mass, &decayInfo.daug1Mass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug2Mass, &decayInfo.daug2Mass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug3Mass, &decayInfo.daug3Mass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_meson_radius, &decayInfo.meson_radius, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    static int cacheCount = 0;
    cacheToUse            = cacheCount++;
    registerConstant(decayInfo.resonances.size());
    registerConstant(cacheToUse);

    for(auto &resonance : decayInfo.resonances) {
        // Not going to use amp_imag, but need a dummy index so the resonance size will be consistent.
        // pindices.push_back(resonance->getFunctionIndex());
        // pindices.push_back(resonance->getParameterIndex());
        // resonance->setConstantIndex(cIndex);
        components.push_back(resonance);
    }

    components.push_back(efficiency);

    initialize();

    redoIntegral = new bool[decayInfo.resonances.size()];
    cachedMasses = new fptype[decayInfo.resonances.size()];
    cachedWidths = new fptype[decayInfo.resonances.size()];
    integrals    = new double[decayInfo.resonances.size()];

    for(int i = 0; i < decayInfo.resonances.size(); ++i) {
        redoIntegral[i] = true;
        cachedMasses[i] = -1;
        cachedWidths[i] = -1;
        integrals[i]    = 0;
    }

    integrators = new SpecialIncoherentIntegrator *[decayInfo.resonances.size()];
    calculators = new SpecialIncoherentResonanceCalculator *[decayInfo.resonances.size()];

    for(int i = 0; i < decayInfo.resonances.size(); ++i) {
        integrators[i] = new SpecialIncoherentIntegrator(parameters, i);
        calculators[i] = new SpecialIncoherentResonanceCalculator(parameters, i);
    }

    addSpecialMask(PdfBase::ForceSeparateNorm);
}

__host__ void IncoherentSumPdf::recursiveSetIndices() {
    GET_FUNCTION_ADDR(ptr_to_incoherent);

    GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName(), "ptr_to_incoherent");
    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx                               = num_device_functions++;

    populateArrays();

    // save our efficiency function.  Resonance's are saved first, then the efficiency function.  Take -1 as efficiency!
    efficiencyFunction = num_device_functions - 1;
}
__host__ void IncoherentSumPdf::setDataSize(unsigned int dataSize, unsigned int evtSize) {
    // Default 3 is m12, m13, evtNum
    totalEventSize = evtSize;
    if(totalEventSize < 3)
        throw GooFit::GeneralError("totalEventSize {} must be 3 or more", totalEventSize);

    if(cachedResonances) {
        delete cachedResonances;
    }

    numEntries       = dataSize;
    cachedResonances = new thrust::device_vector<fpcomplex>(dataSize * decayInfo.resonances.size());
    void *dummy      = thrust::raw_pointer_cast(cachedResonances->data());
    MEMCPY_TO_SYMBOL(
        cResonanceValues, &dummy, sizeof(fpcomplex *), cacheToUse * sizeof(fpcomplex *), cudaMemcpyHostToDevice);
    setForceIntegrals();
}

__host__ fptype IncoherentSumPdf::normalize() const {
    recursiveSetNormalisation(1); // Not going to normalize efficiency,
    // so set normalisation factor to 1 so it doesn't get multiplied by zero.
    // Copy at this time to ensure that the SpecialCalculators, which need the efficiency,
    // don't get zeroes through multiplying by the normFactor.
    MEMCPY_TO_SYMBOL(
        d_normalisations, host_normalisations, totalNormalisations * sizeof(fptype), 0, cudaMemcpyHostToDevice);

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

    // Check if efficiency changes force redoing the integrals.
    if(efficiency->parametersChanged()) {
        forceRedoIntegrals = true;
    }

    // Check for changed masses or forced integral redo.
    for(unsigned int i = 0; i < decayInfo.resonances.size(); ++i) {
        redoIntegral[i] = forceRedoIntegrals;

        if(!(decayInfo.resonances[i]->parametersChanged()))
            continue;

        redoIntegral[i] = true;
    }

    forceRedoIntegrals = false;

    thrust::constant_iterator<fptype *> arrayAddress(dalitzNormRange);
    thrust::counting_iterator<int> binIndex(0);

    // NB, SpecialIncoherentResonanceCalculator assumes that fit is unbinned!
    // And it needs to know the total event size, not just observables
    // for this particular PDF component.
    thrust::constant_iterator<fptype *> dataArray(dev_event_array);
    thrust::constant_iterator<int> eventSize(totalEventSize);
    thrust::counting_iterator<int> eventIndex(0);

    for(int i = 0; i < decayInfo.resonances.size(); ++i) {
        if(redoIntegral[i]) {
            calculators[i]->setIncoherentIndex(getFunctionIndex());
            calculators[i]->setResonanceIndex(decayInfo.resonances[i]->getFunctionIndex());
            thrust::transform(
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, arrayAddress, eventSize)),
                strided_range<thrust::device_vector<fpcomplex>::iterator>(
                    cachedResonances->begin() + i, cachedResonances->end(), decayInfo.resonances.size())
                    .begin(),
                *(calculators[i]));

            integrators[i]->setIncoherentIndex(getFunctionIndex());
            integrators[i]->setEfficiencyIndex(efficiencyFunction);
            integrators[i]->setResonanceIndex(decayInfo.resonances[i]->getFunctionIndex());
            fptype dummy = 0;
            static thrust::plus<fptype> cudaPlus;
            integrals[i] = thrust::transform_reduce(
                thrust::make_zip_iterator(thrust::make_tuple(binIndex, arrayAddress)),
                thrust::make_zip_iterator(thrust::make_tuple(binIndex + totalBins, arrayAddress)),
                *(integrators[i]),
                dummy,
                cudaPlus);
        }
    }

    // End of time-consuming integrals and caching of BWs over Dalitz plot.

    fptype ret = 0;

    for(unsigned int i = 0; i < decayInfo.resonances.size(); ++i) {
        // int param_i      = parameters + resonanceOffset_incoherent + resonanceSize * i;
        fptype amplitude = host_parameters[parametersIdx + i + 1];
        ret += amplitude * integrals[i];
    }

    double binSizeFactor = 1;
    binSizeFactor *= _m12.getBinSize();
    binSizeFactor *= _m13.getBinSize();
    ret *= binSizeFactor;

    host_normalisations[normalIdx + 1] = 1.0 / ret;
    return ret;
}

SpecialIncoherentIntegrator::SpecialIncoherentIntegrator(int pIdx, unsigned int ri)
    : resonance_i(ri)
    , parameters(pIdx) {}

__device__ fptype SpecialIncoherentIntegrator::operator()(thrust::tuple<int, fptype *> t) const {
    // Returns integral of specific BW over Dalitz plot, to be cached and
    // multiplied by rapidly-changing amplitude.

    // Bin index, base address [lower, upper,getNumBins]
    // Notice that this is basically MetricTaker::operator (binned) with the special-case knowledge
    // that event size is two, and that the function to call is getResonanceAmplitude.

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

    while(pc.funcIdx < incoherentSum)
        pc.incrementIndex();

    if(!inDalitz(binCenterM12, binCenterM13, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass))
        return 0;

    int id_m12 = RO_CACHE(pc.observables[pc.observableIdx + 1]);
    int id_m13 = RO_CACHE(pc.observables[pc.observableIdx + 2]);

    int num_res = RO_CACHE(pc.constants[pc.constantIdx + 5]);

    // int parameter_i
    //    = parIndexFromResIndex_incoherent(resonance_i); // Find position of this resonance relative to TDDP start
    // unsigned int functn_i = indices[parameter_i + 2];
    // unsigned int params_i = indices[parameter_i + 3];
    fptype m23 = c_motherMass * c_motherMass + c_daug1Mass * c_daug1Mass + c_daug2Mass * c_daug2Mass
                 + c_daug3Mass * c_daug3Mass - binCenterM12 - binCenterM13;

    while(pc.funcIdx < num_res)
        pc.incrementIndex();

    fpcomplex ret = getResonanceAmplitude(binCenterM12, binCenterM13, m23, pc);

    while(pc.funcIdx < efficiency)
        pc.incrementIndex();

    // unsigned int numResonances = indices[2];
    fptype fakeEvt[10]; // Need room for many observables in case m12 or m13 were assigned a high index in an
                        // event-weighted fit.
    fakeEvt[0]      = 2;
    fakeEvt[id_m12] = binCenterM12;
    fakeEvt[id_m13] = binCenterM13;
    // int effFunctionIdx                   = parIndexFromResIndex_incoherent(numResonances);
    fptype eff = callFunction(fakeEvt, pc);

    return thrust::norm(ret) * eff;
}

SpecialIncoherentResonanceCalculator::SpecialIncoherentResonanceCalculator(int pIdx, unsigned int res_idx)
    : resonance_i(res_idx)
    , parameters(pIdx) {}

__device__ fpcomplex SpecialIncoherentResonanceCalculator::operator()(thrust::tuple<int, fptype *, int> t) const {
    // Returns the BW, or other resonance function, for a specific resonance.
    // Is special because the value is expected to change slowly, so it's
    // useful to cache the result.
    int evtNum  = thrust::get<0>(t);
    fptype *evt = thrust::get<1>(t) + (evtNum * thrust::get<2>(t));

    // unsigned int *indices = paramIndices + parameters; // Jump to TDDP position within parameters array
    ParameterContainer pc;

    while(pc.funcIdx < incoherentSum)
        pc.incrementIndex();

    int id_m12 = RO_CACHE(pc.observables[pc.observableIdx + 1]);
    int id_m13 = RO_CACHE(pc.observables[pc.observableIdx + 2]);

    fptype m12 = evt[id_m12];
    fptype m13 = evt[id_m13];

    if(!inDalitz(m12, m13, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass))
        return {0.0, 0.0};

    fptype m23 = c_motherMass * c_motherMass + c_daug1Mass * c_daug1Mass + c_daug2Mass * c_daug2Mass
                 + c_daug3Mass * c_daug3Mass - m12 - m13;

    while(pc.funcIdx < resonance_i)
        pc.incrementIndex();

    // int parameter_i
    //    = parIndexFromResIndex_incoherent(resonance_i); // Find position of this resonance relative to TDDP start
    // unsigned int functn_i       = indices[parameter_i + 2];
    // unsigned int params_i       = indices[parameter_i + 3];
    fpcomplex ret = getResonanceAmplitude(m12, m13, m23, pc);

    return ret;
}

} // namespace GooFit
