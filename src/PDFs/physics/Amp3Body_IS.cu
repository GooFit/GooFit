#include <goofit/Error.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/Amp3Body_IS.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/PDFs/physics/detail/SpecialIncoherentIntegrator.h>
#include <goofit/PDFs/physics/detail/SpecialIncoherentResonanceCalculator.h>

#include <thrust/complex.h>
#include <thrust/transform_reduce.h>

namespace GooFit {

// Offset of the first resonance into the parameter index array.
// Notice that this is different from the Amp3Body_TD case because there's no time information.
// In particular the offset consists of nP, constant index, number of resonances, and cache index.
const int resonanceOffset_incoherent = 4;

__device__ fpcomplex *cResonanceValues[10];

__device__ inline auto parIndexFromResIndex_incoherent(int resIndex) -> int {
    return resonanceOffset_incoherent + resIndex * resonanceSize;
}

__device__ auto device_incoherent(fptype *evt, ParameterContainer &pc) -> fptype {
    // Calculates the incoherent sum over the resonances.
    int evtId   = pc.getObservable(2);
    auto evtNum = static_cast<int>(floor(0.5 + RO_CACHE(evt[evtId])));

    fptype ret                 = 0;
    unsigned int numResonances = pc.getConstant(4);
    unsigned int cacheToUse    = pc.getConstant(5);

    for(int i = 0; i < numResonances; ++i) {
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

__host__ Amp3Body_IS::Amp3Body_IS(
    std::string n, Observable m12, Observable m13, EventNumber eventNumber, DecayInfo3 decay, GooPdf *eff)
    : Amp3BodyBase("Amp3Pdf_IS", n, m12, m13, eventNumber)
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

    setSeparateNorm();

    registerFunction("ptr_to_incoherent", ptr_to_incoherent);
}

__host__ void Amp3Body_IS::populateArrays() {
    PdfBase::populateArrays();

    // save our efficiency function.  Resonance's are saved first, then the efficiency function.  Take -1 as efficiency!
    efficiencyFunction = host_function_table.size() - 1;
}
__host__ void Amp3Body_IS::setDataSize(unsigned int dataSize, unsigned int evtSize) {
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

__host__ auto Amp3Body_IS::normalize() -> fptype {
    recursiveSetNormalization(1.0); // Not going to normalize efficiency,
    // so set normalization factor to 1 so it doesn't get multiplied by zero.
    // Copy at this time to ensure that the SpecialCalculators, which need the efficiency,
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

    host_normalizations[normalIdx + 1] = 1.0 / ret;
    cachedNormalization                = 1.0 / ret;
    return ret;
}

} // namespace GooFit
