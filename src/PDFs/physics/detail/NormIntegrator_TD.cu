#include <goofit/PDFs/physics/detail/NormIntegrator_TD.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/Amp4BodyGlobals.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/PDFs/physics/detail/Dim5.h>
#include <goofit/detail/Complex.h>

#include <thrust/functional.h>

namespace GooFit {

NormIntegrator_TD::NormIntegrator_TD() = default;

__device__ auto NormIntegrator_TD::operator()(thrust::tuple<int, int, fptype *, fpcomplex *> t) const
    -> thrust::tuple<fptype, fptype, fptype, fptype> {
    // unsigned int *indices = paramIndices + _parameters;
    // unsigned int totalAMP = indices[5];

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

    auto AmpAB = AmpA * conj(AmpB);
    return {thrust::norm(AmpA), thrust::norm(AmpB), AmpAB.real(), AmpAB.imag()};
}

} // namespace GooFit
