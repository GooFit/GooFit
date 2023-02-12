#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/detail/SpecialSqDpResonanceIntegrator.h>
#include <goofit/PDFs/physics/Amp3BodySqDP.h>

namespace GooFit {

__device__ auto device_SqDalitzPlot_calcIntegrals(fptype mprime, fptype thetaprime, int res_i, int res_j, ParameterContainer &pc)
    -> fpcomplex {
    // Calculates BW_i(mprime, thetaprime) * BW_j^*(mprime, thetaprime).
    // This calculation is in a separate function so
    // it can be cached. Note that this function expects
    // to be called on a normalization grid, not on
    // observed points, that's why it doesn't use
    // cResonances. No need to cache the values at individual
    // grid points - we only care about totals.
    fptype m12 = calc_m12(mprime,c_motherMass,c_daug1Mass,c_daug2Mass,c_daug3Mass);
    fptype m13 = calc_m13(m12,cos(thetaprime*M_PI), c_motherMass,c_daug1Mass,c_daug2Mass,c_daug3Mass);
    fptype s12 = m12*m12;
    fptype s13 = m13*m13;
    fptype s23 = c_motherMass * c_motherMass + c_daug1Mass * c_daug1Mass + c_daug2Mass * c_daug2Mass
                 + c_daug3Mass * c_daug3Mass - s12 - s13;

    fpcomplex ret;

    if(!inSqDalitz(mprime, thetaprime))
        return ret;

    ParameterContainer ipc = pc;
    while(ipc.funcIdx < res_i)
        ipc.incrementIndex();

    ret = getResonanceAmplitude(s13, s23 , s12, ipc);

    ParameterContainer jpc = pc;
    while(jpc.funcIdx < res_j)
        jpc.incrementIndex();

    ret *= conj(getResonanceAmplitude( s13, s23 , s12, jpc));

    return ret;
}

SpecialSqDpResonanceIntegrator::SpecialSqDpResonanceIntegrator(int pIdx, unsigned int ri, unsigned int rj)
    : resonance_i(ri)
    , resonance_j(rj)
    , parameters(pIdx) {}

__device__ auto SpecialSqDpResonanceIntegrator::operator()(thrust::tuple<int, fptype *, int> t) const -> fpcomplex {
    //(brad): new indexing plan: bin number, function id, parameter id (not required), fptype with actual
    // bins(needed???)
    // Bin index, base address [lower, upper,getNumBins]
    // Notice that this is basically MetricTaker::operator (binned) with the special-case knowledge
    // that event size is two, and that the function to call is dev_DalitzPlot_calcIntegrals.

    int globalBinNumber  = thrust::get<0>(t);
    fptype lowerBoundMPrime = thrust::get<1>(t)[0];
    fptype upperBoundMPrime = thrust::get<1>(t)[1];
    auto numBinsMPrime      = static_cast<int>(floor(thrust::get<1>(t)[2] + 0.5));
    int binNumberMPrime     = globalBinNumber % numBinsMPrime;
    fptype binCenterMPrime  = upperBoundMPrime - lowerBoundMPrime;
    binCenterMPrime /= numBinsMPrime;
    binCenterMPrime *= (binNumberMPrime + 0.5);
    binCenterMPrime += lowerBoundMPrime;

    globalBinNumber /= numBinsMPrime;
    fptype lowerBoundThetaPrime = thrust::get<1>(t)[3];
    fptype upperBoundThetaPrime = thrust::get<1>(t)[4];
    auto numBinsThetaPrime      = static_cast<int>(floor(thrust::get<1>(t)[2] + 0.5));
    fptype binCenterThetaPrime  = upperBoundThetaPrime - lowerBoundThetaPrime;
    binCenterThetaPrime /= numBinsThetaPrime;
    binCenterThetaPrime *= (globalBinNumber + 0.5);
    binCenterThetaPrime += lowerBoundThetaPrime;

    ParameterContainer pc;

    fptype events[10];

    while(pc.funcIdx < dalitz_i)
        pc.incrementIndex();

    fpcomplex ret = device_SqDalitzPlot_calcIntegrals(binCenterMPrime, binCenterThetaPrime, resonance_i, resonance_j, pc);

    // TODO: read id's in in order to set them for the fake event.

    int id_m13 = pc.getObservable(0);
    int id_m23 = pc.getObservable(1);

    // fptype fakeEvt[10]; // Need room for many observables in case mprime or thetaprime were assigned a high index in an
    // event-weighted fit.
    // fakeEvt[0] = 2;
    // fakeEvt[id_m12] = binCenterMPrime;
    // fakeEvt[id_m13] = binCenterThetaPrime;

    fptype m12 = calc_m12(binCenterMPrime,c_motherMass,c_daug1Mass,c_daug2Mass,c_daug3Mass);
    fptype m13 = calc_m13(m12,cos(binCenterThetaPrime*M_PI), c_motherMass,c_daug1Mass,c_daug2Mass,c_daug3Mass);
    fptype s12 = m12*m12;
    fptype s13 = m13*m13;
    fptype s23 = c_motherMass * c_motherMass + c_daug1Mass * c_daug1Mass + c_daug2Mass * c_daug2Mass
                 + c_daug3Mass * c_daug3Mass - s12 - s13;

    events[0]      = 2;
    events[id_m13] = s13;
    events[id_m23] = s23;

    // unsigned int numResonances           = indices[2];
    // int effFunctionIdx                   = parIndexFromResIndex_DP(numResonances);

    // increment until we are on the efficiency function (17)
    int effFunc = thrust::get<2>(t);
    while(pc.funcIdx < effFunc)
        pc.incrementIndex();

    // make sure that the efficiency function does not use RO_CACHE to access events, as this will cause a crash on GPU
    // (accessing stack memory is not allowed with ldg)
    fptype eff = callFunction(events, pc);

    if(m_no_eff)
        eff = 1.;

    // Multiplication by eff, not sqrt(eff), is correct:
    // These complex numbers will not be squared when they
    // go into the integrals. They've been squared already,
    // as it were.
    fptype jacobian = calc_SqDp_Jacobian(binCenterMPrime, binCenterThetaPrime, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass);
    ret*=eff*jacobian;
    // printf("ret %f %f %f %f %f\n",binCenterMPrime, binCenterThetaPrime, ret.real, ret.imag, eff );
    return ret;
}

} // namespace GooFit
