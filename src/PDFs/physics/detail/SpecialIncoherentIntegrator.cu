#include <goofit/PDFs/physics/detail/SpecialIncoherentIntegrator.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/detail/Complex.h>

namespace GooFit {

SpecialIncoherentIntegrator::SpecialIncoherentIntegrator(int pIdx, unsigned int ri)
    : resonance_i(ri)
    , parameters(pIdx) {}

__device__ auto SpecialIncoherentIntegrator::operator()(thrust::tuple<int, fptype *> t) const -> fptype {
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

    int id_m12 = pc.getObservable(0);
    int id_m13 = pc.getObservable(1);

    int num_res = pc.getConstant(4);

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
    auto *fakeEvt = new fptype[10];

    // Need room for many observables in case m12 or m13 were assigned a high index in an
    // event-weighted fit.
    fakeEvt[0]      = 2;
    fakeEvt[id_m12] = binCenterM12;
    fakeEvt[id_m13] = binCenterM13;
    // int effFunctionIdx                   = parIndexFromResIndex_incoherent(numResonances);
    fptype eff = callFunction(fakeEvt, pc);

    delete[] fakeEvt;

    return thrust::norm(ret) * eff;
}

} // namespace GooFit
