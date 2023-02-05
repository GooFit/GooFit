#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/Amp3Body_TD.h>
#include <goofit/PDFs/physics/detail/SpecialDalitzIntegrator.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

#include <thrust/transform_reduce.h>

namespace GooFit {

__device__ auto device_Tddp_calcIntegrals(fptype m12, fptype m13, int res_i, int res_j, ParameterContainer &pc)
    -> ThreeComplex {
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

SpecialDalitzIntegrator::SpecialDalitzIntegrator(int pIdx, unsigned int ri, unsigned int rj)
    : resonance_i(ri)
    , resonance_j(rj)
    , parameters(pIdx) {}

__device__ auto SpecialDalitzIntegrator::operator()(thrust::tuple<int, fptype *, int> t) const -> ThreeComplex {
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

    // TODO: should we pass in the evtSize?
    auto *events = new fptype[10];

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
    // fptype eff = (*(reinterpret_cast<device_function_ptr>(d_function_table[indices[effFunctionIdx]])))(fakeEvt,
    // cudaArray, paramIndices + indices[effFunctionIdx + 1]);
    while(pc.funcIdx < thrust::get<2>(t))
        pc.incrementIndex();

    fptype eff = callFunction(events, pc);

    delete[] events;
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

} // namespace GooFit
