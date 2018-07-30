#include <goofit/PDFs/physics/detail/AmpCalc_TD.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/Amp4BodyGlobals.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/PDFs/physics/detail/EvalVar.h>

namespace GooFit {

AmpCalc_TD::AmpCalc_TD(unsigned int nPerm, unsigned int ampIdx)
    : _nPerm(nPerm)
    , _AmpIdx(ampIdx) {}

__device__ fpcomplex AmpCalc_TD::operator()(thrust::tuple<int, fptype *, int> t) const {
    ParameterContainer pc;

    while(pc.funcIdx < dalitzFuncId)
        pc.incrementIndex();

    unsigned int cacheToUse = pc.getConstant(5);
    unsigned int totalLS    = pc.getConstant(6);
    unsigned int totalSF    = pc.getConstant(7);
    unsigned int totalAMP   = pc.getConstant(8);
    unsigned int offset     = totalLS + totalSF;
    unsigned int numLS      = AmpIndices[totalAMP + _AmpIdx];
    unsigned int numSF      = AmpIndices[totalAMP + _AmpIdx + 1];
    unsigned int evtNum     = thrust::get<0>(t);

    fpcomplex returnVal(0, 0);
    unsigned int SF_step = numSF / _nPerm;
    unsigned int LS_step = numLS / _nPerm;

    for(int i = 0; i < _nPerm; ++i) {
        fpcomplex ret(1, 0);

        for(int j = i * LS_step; j < (i + 1) * LS_step; ++j) {
            int idx = AmpIndices[totalAMP + _AmpIdx + 4 + j];
            ret *= (cResSF_TD[cacheToUse][evtNum * offset + idx]);
            // printf("Lineshape %i = (%.7g, %.7g)\n", j, (cResSF_TD[cacheToUse][evtNum*offset + AmpIndices[totalAMP +
            // _AmpIdx + 4 + j]]).real, (cResSF_TD[cacheToUse][evtNum*offset + AmpIndices[totalAMP + _AmpIdx + 4 +
            // j]]).imag);
        }

        // printf("Lineshape Product = (%.7g, %.7g)\n", ret.real, ret.imag);
        for(int j = i * SF_step; j < (i + 1) * SF_step; ++j) {
            int idx = AmpIndices[totalAMP + _AmpIdx + 4 + numLS + j];
            ret *= (cResSF_TD[cacheToUse][evtNum * offset + totalLS + idx].real());
            // printf(" SF = %.7g\n", (cResSF_TD[cacheToUse][evtNum*offset + totalLS + AmpIndices[totalAMP + _AmpIdx + 4
            // + numLS + j]].real));
        }

        // printf("Lineshape Product * SF = (%.7g, %.7g)\n", ret.real, ret.imag);

        returnVal += ret;
    }

    returnVal *= (1 / sqrt(static_cast<fptype>(_nPerm)));
    // printf("Amplitude Value = (%.7g, %.7g)\n", returnVal.real(), returnVal.imag());
    return returnVal;
}

} // namespace GooFit
