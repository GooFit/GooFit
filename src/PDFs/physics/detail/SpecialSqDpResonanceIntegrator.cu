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

//    printf("pc = %d \n",pc.funcIdx); == 1

    fpcomplex ret(0.,0.);

    if(!inSqDalitz(mprime, thetaprime))
        return ret;


    fptype m12 = calc_m12(mprime,c_motherMass,c_daug1Mass,c_daug2Mass,c_daug3Mass);
    fptype m13 = calc_m13(m12,cos(thetaprime*M_PI), c_motherMass,c_daug1Mass,c_daug2Mass,c_daug3Mass);
    fptype s12 = m12*m12;
    fptype s13 = m13*m13;
    fptype s23 = c_motherMass * c_motherMass + c_daug1Mass * c_daug1Mass + c_daug2Mass * c_daug2Mass
                 + c_daug3Mass * c_daug3Mass - s12 - s13;

    

    

    if(!inDalitz2(s13, s23,c_motherMass,c_daug1Mass,c_daug2Mass,c_daug3Mass ))
        return ret;

    ParameterContainer ipc = pc;
    while(ipc.funcIdx < res_i)
        ipc.incrementIndex();

    ret = getResonanceAmplitude(s13, s23 , s12 , ipc);

    

    ParameterContainer jpc = pc;
    while(jpc.funcIdx < res_j)
        jpc.incrementIndex();
    

    ret *= conj(getResonanceAmplitude(s13, s23 , s12 , jpc));

    return ret;
}

SpecialSqDpResonanceIntegrator::SpecialSqDpResonanceIntegrator(int pIdx, unsigned int ri, unsigned int rj)
    : resonance_i(ri)
    , resonance_j(rj)
    , parameters(pIdx) {}

__device__ auto SpecialSqDpResonanceIntegrator::operator()(thrust::tuple<int, fptype *, int, int> t) const -> fpcomplex {
   
    int evtNum  = thrust::get<0>(t);
    fptype *evt = thrust::get<1>(t) + (evtNum * thrust::get<2>(t));

    ParameterContainer pc;

    fptype events[20];

    while(pc.funcIdx < dalitz_i)
        pc.incrementIndex();

    int id_mprime = pc.getObservable(0);
    int id_thetaprime = pc.getObservable(1);

    fptype mprime = RO_CACHE(evt[id_mprime]);
    fptype thetaprime = RO_CACHE(evt[id_thetaprime]);

    if(thetaprime>0.5)
        thetaprime = 1.0-thetaprime;

    if(!inSqDalitz(mprime, thetaprime))
        return fpcomplex(0.,0.);

 
    fpcomplex ret = device_SqDalitzPlot_calcIntegrals(mprime, thetaprime, resonance_i, resonance_j, pc);

    // fptype m12 = calc_m12(mprime,c_motherMass,c_daug1Mass,c_daug2Mass,c_daug3Mass);
    // fptype m13 = calc_m13(m12,cos(thetaprime*M_PI), c_motherMass,c_daug1Mass,c_daug2Mass,c_daug3Mass);
    // fptype s12 = m12*m12;
    // fptype s13 = m13*m13;
    // fptype s23 = c_motherMass * c_motherMass + c_daug1Mass * c_daug1Mass + c_daug2Mass * c_daug2Mass
    //              + c_daug3Mass * c_daug3Mass - s12 - s13;


    events[0] = 2;
    events[id_mprime] = mprime;
    events[id_thetaprime] = thetaprime;

    int effFunc = thrust::get<3>(t);

    while(pc.funcIdx < effFunc)
         pc.incrementIndex();

    fptype eff = callFunction(events, pc);

    // printf("mp=%f  th=%f oi pj=%d  eff=%f \n",mprime,thetaprime,pc.funcIdx, eff);
    

    // Multiplication by eff, not sqrt(eff), is correct:
    // These complex numbers will not be squared when they
    // go into the integrals. They've been squared already,
    // as it were.
    fptype jacobian = calc_SqDp_Jacobian(mprime, thetaprime, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass);
    if(m_no_eff)
        return ret*jacobian;
    else
        return ret*eff;

   
    // printf("ret %f %f %f %f %f\n",binCenterMPrime, binCenterThetaPrime, ret.real, ret.imag, eff );
    // return ret;
}

} // namespace GooFit
