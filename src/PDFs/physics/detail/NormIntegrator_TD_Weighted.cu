#include <goofit/PDFs/physics/detail/NormIntegrator_TD_Weighted.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/Amp4BodyGlobals.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/PDFs/physics/detail/Dim5.h>
#include <goofit/detail/Complex.h>
#include <goofit/PDFs/physics/MixingTimeResolution.h>

#include <thrust/functional.h>

namespace GooFit {

NormIntegrator_TD_Weighted::NormIntegrator_TD_Weighted() = default;

// TODO: This entire function essentially just evaluates Amp4Body_TD. We
// should really be able to avoid this.
__device__ fptype NormIntegrator_TD_Weighted::operator()(
    thrust::tuple<int, int, fptype *, fptype *, fptype *, fptype *, fpcomplex *> t) const {
    ParameterContainer pc;

    while(pc.funcIdx < dalitzFuncId)
        pc.incrementIndex();

    unsigned int totalAMP     = pc.getConstant(8);
    unsigned int totalSFandLS = pc.getConstant(10);

    unsigned int evtNum   = thrust::get<0>(t);
    unsigned int MCevents = thrust::get<1>(t);
    fptype dtime          = thrust::get<2>(t)[evtNum];
    fptype sigma          = thrust::get<3>(t)[evtNum];
    fptype weight         = thrust::get<4>(t)[evtNum];
    fptype *SFnorm        = thrust::get<5>(t) + evtNum;
    fpcomplex *LSnorm     = thrust::get<6>(t) + evtNum;

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

    auto AmpAB = AmpA * conj(AmpB);

    fptype tau     = pc.getParameter(0);
    fptype xmixing = pc.getParameter(1);
    fptype ymixing = pc.getParameter(2);
    fptype term1   = thrust::norm(AmpA) + thrust::norm(AmpB);
    fptype term2   = thrust::norm(AmpA) - thrust::norm(AmpB);

    // Find the resolution function.
    pc.incrementIndex();
    for(int i = 0; i < totalSFandLS; i++)
        pc.incrementIndex();

    fptype ret = (*(reinterpret_cast<device_resfunction_ptr>(d_function_table[pc.funcIdx])))(
        term1, term2, AmpAB.real(), AmpAB.imag(), tau, dtime, xmixing, ymixing, sigma, pc);

    ret *= weight;

    return ret;
}

} // namespace GooFit
