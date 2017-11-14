#include "goofit/PDFs/basic/KinLimitBWPdf.h"
#include "goofit/Variable.h"

namespace GooFit {

__device__ fptype getMomentum(const fptype &mass, const fptype &pimass, const fptype &d0mass) {
    if(mass <= 0)
        return 0;

    double lambda = mass * mass - pimass * pimass - d0mass * d0mass;
    lambda *= lambda;
    lambda -= 4 * pimass * pimass * d0mass * d0mass;

    if(lambda <= 0)
        return 0;

    return sqrt(0.5 * lambda / mass);
}

__device__ fptype bwFactor(const fptype &momentum) {
    // 2.56 = 1.6^2, comes from radius for spin-1 particle
    return 1 / sqrt(1.0 + 2.56 * momentum * momentum);
}

__device__ fptype device_KinLimitBW(fptype *evt, ParameterContainer &pc) {
    int id = RO_CACHE(pc.observables[pc.observableIdx + 1]);

    fptype x      = evt[id];
    fptype mean   = RO_CACHE(pc.parameters[pc.parameterIdx + 1]);
    fptype width  = RO_CACHE(pc.parameters[pc.parameterIdx + 2]);
    fptype d0mass = RO_CACHE(pc.constants[pc.constantIdx + 1]);
    fptype pimass = RO_CACHE(pc.constants[pc.constantIdx + 2]);

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
    : GooPdf(n, _x) {
    registerParameter(mean);
    registerParameter(width);

    registerConstant(1.8645);
    registerConstant(0.13957);

    initialize();
}

__host__ void KinLimitBWPdf::recursiveSetIndices() {
    GET_FUNCTION_ADDR(ptr_to_KinLimitBW);

    GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName(), "ptr_to_KinLimitBW");
    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx                               = num_device_functions++;

    populateArrays();
}

__host__ void KinLimitBWPdf::setMasses(fptype bigM, fptype smallM) {
    constantsList[0] = bigM;
    constantsList[1] = smallM;
    // MEMCPY_TO_SYMBOL(functorConstants, constants, 2 * sizeof(fptype), cIndex * sizeof(fptype),
    // cudaMemcpyHostToDevice);
}
} // namespace GooFit
