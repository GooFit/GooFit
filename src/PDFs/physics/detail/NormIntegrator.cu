#include <goofit/PDFs/physics/detail/NormIntegrator.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/Amp4BodyGlobals.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/PDFs/physics/detail/Dim5.h>
#include <goofit/detail/Complex.h>

#include <thrust/functional.h>

namespace GooFit {

NormIntegrator::NormIntegrator() = default;

__device__ auto NormIntegrator::operator()(thrust::tuple<int, int, fptype *, fpcomplex *> t) const -> fptype {
    ParameterContainer pc;

    while(pc.funcIdx < dalitzFuncId)
        pc.incrementIndex();

    unsigned int totalAMP = pc.getConstant(8);

    unsigned int evtNum   = thrust::get<0>(t);
    unsigned int MCevents = thrust::get<1>(t);
    fptype *SFnorm        = thrust::get<2>(t) + evtNum;
    fpcomplex *LSnorm     = thrust::get<3>(t) + evtNum;

    fpcomplex returnVal(0, 0);

    for(int amp = 0; amp < totalAMP; ++amp) {
        unsigned int ampidx  = AmpIndices[amp];
        unsigned int numLS   = AmpIndices[totalAMP + ampidx];
        unsigned int numSF   = AmpIndices[totalAMP + ampidx + 1];
        unsigned int nPerm   = AmpIndices[totalAMP + ampidx + 2];
        unsigned int SF_step = numSF / nPerm;
        unsigned int LS_step = numLS / nPerm;
        fpcomplex ret2(0, 0);
        // printf("%i, %i, %i, %i, %i, %i, %i, %i, %i, %f\n",ampidx, amp, numLS, numSF, nPerm,AmpIndices[totalAMP +
        // ampidx + 3 + 0], AmpIndices[totalAMP + ampidx + 3 + 1], AmpIndices[totalAMP + ampidx + 3 + 2],
        // AmpIndices[totalAMP + ampidx + 3 + 3], (1/sqrt((fptype)(nPerm))) );

        for(int j = 0; j < nPerm; ++j) {
            fpcomplex ret(1, 0);

            for(int i = j * LS_step; i < (j + 1) * LS_step; ++i) {
                fpcomplex matrixelement(LSnorm[AmpIndices[totalAMP + ampidx + 3 + i] * MCevents]);
                // printf("Norm BW %i, %.5g, %.5g\n",AmpIndices[totalAMP + ampidx + 3 + i] , matrixelement.real,
                // matrixelement.imag);
                ret *= matrixelement;
            }

            for(int i = j * SF_step; i < (j + 1) * SF_step; ++i) {
                fptype matrixelement = (SFnorm[AmpIndices[totalAMP + ampidx + 3 + numLS + i] * MCevents]);
                // printf("Norm SF %i, %.5g\n",AmpIndices[totalAMP + ampidx + 3 + i] , matrixelement);
                ret *= matrixelement;
            }

            ret2 += ret;
        }

        fpcomplex amp_C{pc.getParameter(2 * amp + 0), pc.getParameter(2 * amp + 1)};
        ret2 *= (1 / sqrt(static_cast<fptype>(nPerm)));
        // printf("Result Amplitude %i, %.5g, %.5g\n",amp, ret2.real, ret2.imag);
        returnVal += ret2 * amp_C;
    }

    return thrust::norm(returnVal);
}

} // namespace GooFit
