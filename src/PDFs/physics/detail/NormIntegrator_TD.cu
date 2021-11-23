#include <goofit/PDFs/physics/detail/NormIntegrator_TD.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/Amp4BodyGlobals.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/PDFs/physics/detail/Dim5.h>
#include <goofit/detail/Complex.h>
#include <goofit/PDFs/physics/MixingTimeResolution.h>
#include <thrust/functional.h>

namespace GooFit {

 //NormIntegrator_TD::NormIntegrator_TD() = default;
 NormIntegrator_TD::NormIntegrator_TD(bool SpecInt): _SpecInt(SpecInt){}
 __device__ thrust::tuple<fptype, fptype, fptype, fptype> NormIntegrator_TD::
 operator()(thrust::tuple<int, int, fptype *, fpcomplex *,
        mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t> t) const {
    // unsigned int *indices = paramIndices + _parameters;
    // unsigned int totalAMP = indices[5];
    //printf("Norm integrator operator()\n");
    ParameterContainer pc;

    while(pc.funcIdx < dalitzFuncId)
        pc.incrementIndex();

    unsigned int totalAMP = pc.getConstant(8);

    unsigned int evtNum   = thrust::get<0>(t);
    unsigned int MCevents = thrust::get<1>(t);
    fptype *SFnorm        = thrust::get<2>(t) + evtNum;
    fpcomplex *LSnorm     = thrust::get<3>(t) + evtNum;

    fpcomplex AmpA(0, 0);
    fpcomplex AmpB(0, 0);
    fpcomplex amp_A, amp_B;

    int k = 0;

    for(int amp = 0; amp < totalAMP; ++amp) {
        unsigned int ampidx  = AmpIndices[amp];
        unsigned int numLS   = AmpIndices[totalAMP + ampidx];
        unsigned int numSF   = AmpIndices[totalAMP + ampidx + 1];
        unsigned int nPerm   = AmpIndices[totalAMP + ampidx + 2];
        unsigned int flag    = AmpIndices[totalAMP + ampidx + 3];
        unsigned int SF_step = numSF / nPerm;
        unsigned int LS_step = numLS / nPerm;
        fpcomplex ret2(0, 0);
        // printf("%i, %i, %i, %i, %i, %i, %i, %i, %i, %f\n",ampidx, amp, numLS, numSF, nPerm,AmpIndices[totalAMP +
        // ampidx + 4 + 0], AmpIndices[totalAMP + ampidx + 4 + 1], AmpIndices[totalAMP + ampidx + 4 + 2],
        // AmpIndices[totalAMP + ampidx + 4 + 3], (1/sqrt((fptype)(nPerm))) );

        for(int j = 0; j < nPerm; ++j) {
            fpcomplex ret(1, 0);

            for(int i = j * LS_step; i < (j + 1) * LS_step; ++i) {
                fpcomplex matrixelement(LSnorm[AmpIndices[totalAMP + ampidx + 4 + i] * MCevents]);
                // printf("Norm BW %i, %.5g, %.5g\n",AmpIndices[totalAMP + ampidx + 4 + i] , matrixelement.real,
                // matrixelement.imag);
                ret *= matrixelement;
            }

            for(int i = j * SF_step; i < (j + 1) * SF_step; ++i) {
                fptype matrixelement = (SFnorm[AmpIndices[totalAMP + ampidx + 4 + numLS + i] * MCevents]);
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

    fptype _SqWStoRSrate = pc.getParameter(3);
    AmpA *= _SqWStoRSrate;
    //printf("Before if statement\n");
    auto AmpAB = AmpA * conj(AmpB);
    if(_SpecInt){
        //printf("Inside if statement\n");
        fptype _tau          = pc.getParameter(0);
        fptype _xmixing      = pc.getParameter(1);
        fptype _ymixing      = pc.getParameter(2);
        fptype _time = thrust::get<4>(t);
        fptype _eff = thrust::get<5>(t);
        //importance weight 
        fptype _importance_weight = thrust::get<6>(t);
        fptype term1 = thrust::norm(AmpA) + thrust::norm(AmpB);
        fptype term2 = thrust::norm(AmpA) - thrust::norm(AmpB);
  
        unsigned int totalSF_LS = pc.getConstant(10);
        fptype ret = 0.;
      
        //printf("Calling resolution function\n");
        //ret = (*(reinterpret_cast<device_resfunction_ptr>(d_function_table[pc.funcIdx])))(
        //term1, term2, AmpAB.real(), AmpAB.imag(), _tau, _time, _xmixing, _ymixing,0.0, pc);
        /*
        _time /= _tau;
        ret += term1 * cosh(_ymixing * _time);
        ret += term2 * cos(_xmixing * _time);
        ret -= 2 * AmpAB.real() * sinh(_ymixing * _time);
        ret -= 2 * AmpAB.imag()* sin(_xmixing * _time); // Notice sign difference wrt to Mikhail's code, because I have AB* and he has A*B.                                                                                      
        ret *= exp(-_time);
        */;
        //printf("Called resolution function\n");
        //printf("ret step 1: %f\n",ret);
        ret *= _eff;
        //printf("ret step 2: %f\n",ret);
        ret *= _importance_weight;
        //printf("ret step 3: %f\n",ret);

        /*
        fptype integral1 = _tau / (1 - _ymixing * _ymixing);
        fptype integral2 = _tau / (1 + _xmixing * _xmixing);
        fptype integral3 = _ymixing * integral1;
        fptype integral4 = _xmixing * integral2;
        ret += integral2 * (thrust::norm(AmpA) -  thrust::norm(AmpB));
        ret -= 2 * integral3 * AmpAB.real();
        ret -= 2 * integral4 * AmpAB.imag();
        ret *= _eff;
        ret *= _weight;
        ret /= _importance_weight;
        */
  
        return thrust::tuple<fptype,fptype,fptype,fptype>(ret,thrust::norm(AmpA)/_importance_weight,thrust::norm(AmpB)/_importance_weight,1);
      }
      else {
        return {thrust::norm(AmpA), thrust::norm(AmpB), AmpAB.real(), AmpAB.imag()};
        }
}

} // namespace GooFit
