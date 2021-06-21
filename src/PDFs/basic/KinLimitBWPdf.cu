#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/basic/KinLimitBWPdf.h>
#include <goofit/Variable.h>

namespace GooFit {

__device__ auto getMomentum(const fptype &mass, const fptype &pimass, const fptype &d0mass) -> fptype {
    if(mass <= 0)
        return 0;

    double lambda = mass * mass - pimass * pimass - d0mass * d0mass;
    lambda *= lambda;
    lambda -= 4 * pimass * pimass * d0mass * d0mass;

    if(lambda <= 0)
        return 0;

    return sqrt(0.5 * lambda / mass);
}

__device__ auto bwFactor(const fptype &momentum) -> fptype {
    // 2.56 = 1.6^2, comes from radius for spin-1 particle
    return 1 / sqrt(1.0 + 2.56 * momentum * momentum);
}

__device__ auto device_KinLimitBW(fptype *evt, ParameterContainer &pc) -> fptype {
    int id = pc.getObservable(0);

    fptype x      = RO_CACHE(evt[id]);
    fptype mean   = pc.getParameter(0);
    fptype width  = pc.getParameter(1);
    fptype d0mass = pc.getConstant(0);
    fptype pimass = pc.getConstant(1);

    mean += d0mass;
    x += d0mass;

    fptype pUsingRealMass = getMomentum(mean, pimass, d0mass);

    pc.incrementIndex(1, 2, 2, 1, 1);

    if(0 >= pUsingRealMass)
        return 0;

    mean *= mean;
    fptype pUsingX     = getMomentum(x, pimass, d0mass);
    fptype phspfactor  = POW3(pUsingX / pUsingRealMass) * POW2(bwFactor(pUsingX) / bwFactor(pUsingRealMass));
    fptype phspMassSq  = POW2(mean - x * x);
    fptype phspGammaSq = POW2(width * phspfactor);

    fptype ret = (phspfactor * mean * width * width) / (phspMassSq + mean * phspGammaSq);

    return ret;
}

__device__ device_function_ptr ptr_to_KinLimitBW = device_KinLimitBW;

__host__ KinLimitBWPdf::KinLimitBWPdf(std::string n, Observable _x, Variable mean, Variable width)
    : GooPdf("KinLimitBWPdf", n, _x, mean, width) {
    registerConstant(1.8645);
    registerConstant(0.13957);

    registerFunction("ptr_to_KinLimitBW", ptr_to_KinLimitBW);

    initialize();
}

__host__ void KinLimitBWPdf::setMasses(fptype bigM, fptype smallM) {
    constantsList[0] = bigM;
    constantsList[1] = smallM;
}
} // namespace GooFit
