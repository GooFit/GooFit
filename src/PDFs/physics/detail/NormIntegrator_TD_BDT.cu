#include <goofit/PDFs/physics/detail/NormIntegrator_TD_BDT.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/Amp4BodyGlobals.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/PDFs/physics/detail/Dim5.h>
#include <goofit/detail/Complex.h>
#include <goofit/PDFs/physics/MixingTimeResolution.h>

#include <thrust/functional.h>

namespace GooFit {

//NormIntegrator_TD_BDT::NormIntegrator_TD_BDT() = default;

NormIntegrator_TD_BDT::NormIntegrator_TD_BDT(unsigned int CacheIdx)
    :_CacheIdx(CacheIdx){}

//Operator that performs the integration with acceptance. Takes the same arguments as ordinary operator but also expects iterators to vectors of norm_dtime, norm_eff and norm_weight (importance sampling)
__device__ auto NormIntegrator_TD_BDT::operator()(thrust::tuple<int, int, fptype *, fpcomplex *, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t> t) const
    -> thrust::tuple<fptype, fptype, fptype, fptype> {
    // unsigned int *indices = paramIndices + _parameters;
    // unsigned int totalAMP = indices[5];

    ParameterContainer pc;
    unsigned int cacheToUse = pc.getConstant(5);
    cacheToUse = NormIntegrator_TD_BDT::_CacheIdx;
    while(pc.funcIdx < dalitzFuncId)
        pc.incrementIndex();

    unsigned int totalAMP = pc.getConstant(8);
    //printf("Inside NormIntegrator with %i amps\n",totalAMP);
    unsigned int evtNum   = thrust::get<0>(t);
    unsigned int MCevents = thrust::get<1>(t);
    fptype *SFnorm        = thrust::get<2>(t) + evtNum;
    fpcomplex *LSnorm     = thrust::get<3>(t) + evtNum;

    fpcomplex AmpA(0, 0);
    fpcomplex AmpB(0, 0);
    fpcomplex amp_A, amp_B;

    int k = 0;

    for(int amp = 0; amp < totalAMP; ++amp) {
        unsigned int ampidx  = AmpIndices[cacheToUse][amp];
        unsigned int numLS   = AmpIndices[cacheToUse][totalAMP + ampidx];
        unsigned int numSF   = AmpIndices[cacheToUse][totalAMP + ampidx + 1];
        unsigned int nPerm   = AmpIndices[cacheToUse][totalAMP + ampidx + 2];
        unsigned int flag    = AmpIndices[cacheToUse][totalAMP + ampidx + 3];
        //printf("numSF: %i, nPerm: %i\n",numSF,nPerm);
        unsigned int SF_step = numSF / nPerm;
        unsigned int LS_step = numLS / nPerm;
        fpcomplex ret2(0, 0);
        // printf("%i, %i, %i, %i, %i, %i, %i, %i, %i, %f\n",ampidx, amp, numLS, numSF, nPerm,AmpIndices[totalAMP +
        // ampidx + 4 + 0], AmpIndices[totalAMP + ampidx + 4 + 1], AmpIndices[totalAMP + ampidx + 4 + 2],
        // AmpIndices[totalAMP + ampidx + 4 + 3], (1/sqrt((fptype)(nPerm))) );

        for(int j = 0; j < nPerm; ++j) {
            fpcomplex ret(1, 0);

            for(int i = j * LS_step; i < (j + 1) * LS_step; ++i) {
                fpcomplex matrixelement(LSnorm[AmpIndices[cacheToUse][totalAMP + ampidx + 4 + i] * MCevents]);
                // printf("Norm BW %i, %.5g, %.5g\n",AmpIndices[totalAMP + ampidx + 4 + i] , matrixelement.real,
                // matrixelement.imag);
                ret *= matrixelement;
            }

            for(int i = j * SF_step; i < (j + 1) * SF_step; ++i) {
                fptype matrixelement = (SFnorm[AmpIndices[cacheToUse][totalAMP + ampidx + 4 + numLS + i] * MCevents]);
                // printf("Norm SF %i, %.5g\n",AmpIndices[totalAMP + ampidx + 4 + i] , matrixelement);
                ret *= matrixelement;
            }

            ret2 += ret;
        }

        ret2 *= (1 / sqrt(static_cast<fptype>(nPerm)));
        // printf("Result Amplitude %i, %i, %.5g, %.5g\n",flag, amp, ret2.real, ret2.imag);

        switch(flag) {
        case 0:
            amp_A = fpcomplex(pc.getParameter(4 + 2 * (amp + k)), pc.getParameter(5 + 2 * (amp + k)));
            AmpA += ret2 * amp_A;
            break;

        case 1:
            amp_B = fpcomplex(pc.getParameter(4 + 2 * (amp + k)), pc.getParameter(5 + 2 * (amp + k)));
            AmpB += ret2 * amp_B;
            break;

        case 2:
            amp_A = fpcomplex(pc.getParameter(4 + 2 * (amp + k)), pc.getParameter(5 + 2 * (amp + k)));
            AmpA += ret2 * amp_A;
            ++k;
            amp_B = fpcomplex(pc.getParameter(4 + 2 * (amp + k)), pc.getParameter(5 + 2 * (amp + k)));
            AmpB += ret2 * amp_B;
            break;
        }
    }

    auto tmpA = AmpA; //get value of AmpA before multiplying by SqWStoRS rate
    fptype _SqWStoRSrate = pc.getParameter(3);
    AmpA *= _SqWStoRSrate;

    auto AmpAB = AmpA * conj(AmpB);

    
    //This is the special part that is different

    fptype _tau          = pc.getParameter(0);
    fptype _xmixing      = pc.getParameter(1);
    fptype _ymixing      = pc.getParameter(2);
    fptype _time         = thrust::get<4>(t);
    fptype _eff          = thrust::get<5>(t);
    fptype _weight       = thrust::get<6>(t);

    //printf("_tau: %.7g, _xmixing: %.7g, _ymixing: %.7g\n",_tau,_xmixing,_ymixing);

    fptype term1                  = thrust::norm(AmpA) + thrust::norm(AmpB);
    fptype term2                  = thrust::norm(AmpA) - thrust::norm(AmpB);
    thrust::complex<fptype> term3 = AmpA * thrust::conj(AmpB);

    //increment pc to get correct index for resolution function
    unsigned int totalSF_LS = pc.getConstant(10);
    pc.incrementIndex();
    for(int i = 0; i < totalSF_LS; i++)
        pc.incrementIndex();
    //printf("Calling res function\n");
    fptype ret = (*(reinterpret_cast<device_resfunction_ptr>(d_function_table[pc.funcIdx])))(
        term1, term2, term3.real(), term3.imag(), _tau, _time, _xmixing, _ymixing, 0., pc); //assume no resolution (= 0), for now.
    //printf("Correcting return value");
    ret *= _eff;
    ret /= _weight;
    if(std::isnan(ret)){
      //printf("return value from normintegrator_TD: %.7g\n",ret);    
        ret = 0; //This shouldn't be done but this is a temporary workaround!!!!!!
    }
    
    return thrust::tuple<fptype, fptype, fptype, fptype>(ret, thrust::norm(tmpA)/_weight, thrust::norm(AmpB)/_weight, 1);
}


} // namespace GooFit
