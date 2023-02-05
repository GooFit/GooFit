#include <goofit/PDFs/physics/detail/AmpCalc.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/Amp4BodyGlobals.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/PDFs/physics/detail/Dim5.h>

namespace GooFit {

AmpCalc::AmpCalc(unsigned int nPerm, unsigned int amp)
    : _nPerm(nPerm)
    , _AmpIdx(amp) {}

__device__ auto AmpCalc::operator()(thrust::tuple<int, fptype *, int> t) const -> fpcomplex {
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
        fpcomplex tmp(1, 0);

        for(int j = i * LS_step; j < (i + 1) * LS_step; ++j) {
            int idx = AmpIndices[totalAMP + _AmpIdx + 3 + j];
            tmp     = (cResSF[cacheToUse][evtNum * offset + idx]);
            ret *= tmp;
            // printf("Lineshape = (%.7g, %.7g)\n", tmp.real, tmp.imag);
        }

        // printf("Lineshape Product = (%.7g, %.7g)\n", ret.real, ret.imag);
        for(int j = i * SF_step; j < (i + 1) * SF_step; ++j) {
            int idx = AmpIndices[totalAMP + _AmpIdx + 3 + numLS + j];
            tmp     = (cResSF[cacheToUse][evtNum * offset + totalLS + idx].real());
            ret *= tmp;
            // printf("SF = (%.7g, %.7g)\n", tmp.real, tmp.imag);
        }

        // printf("Lineshape Product * SF = (%.7g, %.7g)\n", ret.real, ret.imag);

        returnVal += ret;
    }

    returnVal *= (1 / sqrt(static_cast<fptype>(_nPerm)));
    // printf("Amplitude Value = (%.7g, %.7g)\n", returnVal.real, returnVal.imag);
    return returnVal;
}

} // namespace GooFit
