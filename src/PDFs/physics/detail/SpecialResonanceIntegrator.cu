#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/detail/SpecialResonanceIntegrator.h>

namespace GooFit {

__device__ auto device_DalitzPlot_calcIntegrals(fptype m13, fptype m23, int res_i, int res_j, ParameterContainer &pc)
    -> fpcomplex {
    // Calculates BW_i(m13, m23) * BW_j^*(m13, m23).
    // This calculation is in a separate function so
    // it can be cached. Note that this function expects
    // to be called on a normalization grid, not on
    // observed points, that's why it doesn't use
    // cResonances. No need to cache the values at individual
    // grid points - we only care about totals.
    fptype motherMass = c_motherMass; // RO_CACHE(pc.constants[pc.constantIdx + 4]);
    fptype daug1Mass  = c_daug1Mass;  // RO_CACHE(pc.constants[pc.constantIdx + 5]);
    fptype daug2Mass  = c_daug2Mass;  // RO_CACHE(pc.constants[pc.constantIdx + 6]);
    fptype daug3Mass  = c_daug3Mass;  // RO_CACHE(pc.constants[pc.constantIdx + 7]);


    fpcomplex ret;

    if(!inDalitz2(m13, m23, motherMass, daug1Mass, daug2Mass, daug3Mass)){
        return fpcomplex(0.0,0.0);
    }
    fptype m12
        = motherMass * motherMass + daug1Mass * daug1Mass + daug2Mass * daug2Mass + daug3Mass * daug3Mass - m13 - m23;

    ParameterContainer ipc = pc;
    while(ipc.funcIdx < res_i)
        ipc.incrementIndex();

    ret = getResonanceAmplitude(m13, m23, m12, ipc);

    ParameterContainer jpc = pc;
    while(jpc.funcIdx < res_j)
        jpc.incrementIndex();

    ret *= conj(getResonanceAmplitude(m13, m23, m12, jpc));

    return ret;
}

SpecialResonanceIntegrator::SpecialResonanceIntegrator(int pIdx, unsigned int ri, unsigned int rj)
    : resonance_i(ri)
    , resonance_j(rj)
    , parameters(pIdx) {}

__device__ auto SpecialResonanceIntegrator::operator()(thrust::tuple<int, fptype *, int> t) const -> thrust::tuple<fpcomplex, fpcomplex> {
    //(brad): new indexing plan: bin number, function id, parameter id (not required), fptype with actual
    // bins(needed???)
    // Bin index, base address [lower, upper,getNumBins]
    // Notice that this is basically MetricTaker::operator (binned) with the special-case knowledge
    // that event size is two, and that the function to call is dev_DalitzPlot_calcIntegrals.

    int globalBinNumber  = thrust::get<0>(t);
    fptype lowerBoundM13 = thrust::get<1>(t)[0];
    fptype upperBoundM13 = thrust::get<1>(t)[1];
    int numBinsM13      = static_cast<int>(floor(thrust::get<1>(t)[2] + 0.5));
    auto binNumberM13     = globalBinNumber % numBinsM13;
    fptype binCenterM13  = upperBoundM13 - lowerBoundM13;
    binCenterM13 /= numBinsM13;
    binCenterM13 *= (binNumberM13 + 0.5);
    binCenterM13 += lowerBoundM13;

    //printf("%d %f %f %d %d %f\n",globalBinNumber,lowerBoundM12,upperBoundM12,numBinsM12,binNumberM12,binCenterM12);

    globalBinNumber /= numBinsM13;
    fptype lowerBoundM23 = thrust::get<1>(t)[3];
    fptype upperBoundM23 = thrust::get<1>(t)[4];
    auto numBinsM23      = static_cast<int>(floor(thrust::get<1>(t)[5] + 0.5));
    fptype binCenterM23  = upperBoundM23 - lowerBoundM23;
    binCenterM23 /= numBinsM23;
    binCenterM23 *= (globalBinNumber + 0.5);
    binCenterM23 += lowerBoundM23;

    ParameterContainer pc;

    fptype events[3];

    while(pc.funcIdx < dalitz_i)
        pc.incrementIndex();

    
    if(c_SymDp && binCenterM13>binCenterM23)
        return thrust::make_tuple(fpcomplex(0.,0.),fpcomplex(0.,0.));

    fpcomplex ret = device_DalitzPlot_calcIntegrals(binCenterM13, binCenterM23, resonance_i, resonance_j, pc);

    // TODO: read id's in in order to set them for the fake event.

    int id_m13 = pc.getObservable(0);
    int id_m23 = pc.getObservable(1);

    // fptype fakeEvt[10]; // Need room for many observables in case m12 or m13 were assigned a high index in an
    // event-weighted fit.
    // fakeEvt[0] = 2;
    // fakeEvt[id_m12] = binCenterM12;
    // fakeEvt[id_m13] = binCenterM13;

    events[0]      = 2;
    events[id_m13] = binCenterM13;
    events[id_m23] = binCenterM23;

    // unsigned int numResonances           = indices[2];
    // int effFunctionIdx                   = parIndexFromResIndex_DP(numResonances);

    // increment until we are on the efficiency function (17)
    int effFunc = thrust::get<2>(t);
    while(pc.funcIdx < effFunc)
        pc.incrementIndex();

    // make sure that the efficiency function does not use RO_CACHE to access events, as this will cause a crash on GPU
    // (accessing stack memory is not allowed with ldg)
    fptype eff = callFunction(events, pc);

    // Multiplication by eff, not sqrt(eff), is correct:
    // These complex numbers will not be squared when they
    // go into the integrals. They've been squared already,
    // as it were.

    return thrust::make_tuple(ret,ret*eff);
   
}

} // namespace GooFit
