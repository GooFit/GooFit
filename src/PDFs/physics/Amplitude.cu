#include <goofit/PDFs/physics/Amplitude.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>

namespace GooFit {

Amplitude::Amplitude(std::string uniqueDecayStr,
                     Variable ar,
                     Variable ai,
                     std::vector<Lineshape *> LS,
                     std::vector<SpinFactor *> SF,
                     unsigned int nPerm)
    : GooPdf(uniqueDecayStr)
    , _uniqueDecayStr(uniqueDecayStr)
    , _ar(ar)
    , _ai(ai)
    , _SF(std::move(SF))
    , _LS(std::move(LS))
    , _nPerm(nPerm) {
    for(auto &lineshape : _LS)
        components.push_back(lineshape);

    for(auto &spinfactor : _SF)
        components.push_back(spinfactor);

    initialize();
}

void Amplitude::recursiveSetIndices() {
    // There isn't a device function, should there be one?
    for(auto &component : components)
        component->recursiveSetIndices();
}

bool Amplitude::operator==(const Amplitude &A) const {
    return _uniqueDecayStr == A._uniqueDecayStr && _ar == A._ar && _ai == A._ai && _LS == A._LS && _SF == A._SF
           && _nPerm == A._nPerm;
}

} // namespace GooFit
