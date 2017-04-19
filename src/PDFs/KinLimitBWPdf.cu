#include "goofit/PDFs/KinLimitBWPdf.h"

EXEC_TARGET fptype getMomentum(const fptype& mass, const fptype& pimass, const fptype& d0mass) {
    if(mass <= 0)
        return 0;

    double lambda = mass*mass - pimass*pimass - d0mass*d0mass;
    lambda *= lambda;
    lambda -= 4*pimass*pimass*d0mass*d0mass;

    if(lambda <= 0)
        return 0;

    return SQRT(0.5*lambda/mass);
}

EXEC_TARGET fptype bwFactor(const fptype& momentum) {
    // 2.56 = 1.6^2, comes from radius for spin-1 particle
    return 1/SQRT(1.0 + 2.56 * momentum*momentum);
}

EXEC_TARGET fptype device_KinLimitBW(fptype* evt, fptype* p, unsigned int* indices) {
    fptype x = evt[RO_CACHE(indices[2 + RO_CACHE(indices[0])])];
    fptype mean  = RO_CACHE(p[RO_CACHE(indices[1])]);
    fptype width = RO_CACHE(p[RO_CACHE(indices[2])]);
    fptype d0mass = RO_CACHE(functorConstants[RO_CACHE(indices[3])+0]);
    fptype pimass = RO_CACHE(functorConstants[RO_CACHE(indices[3])+1]);

    mean += d0mass;
    x += d0mass;

    fptype pUsingRealMass = getMomentum(mean, pimass, d0mass);

    if(0 >= pUsingRealMass)
        return 0;

    mean *= mean;
    fptype pUsingX     = getMomentum(x, pimass, d0mass);
    fptype phspfactor  = POW3(pUsingX / pUsingRealMass) * POW2(bwFactor(pUsingX) / bwFactor(pUsingRealMass));
    fptype phspMassSq  = POW2(mean - x*x);
    fptype phspGammaSq = POW2(width*phspfactor);

    fptype ret = (phspfactor * mean*width*width)/(phspMassSq + mean*phspGammaSq);
#ifdef CUDAPRINT
    /*
    if (((0 == THREADIDX) && (0 == BLOCKIDX) && (callnumber < 10)) || (isnan(ret)))
        cuPrintf("KinLimitBW %f %f %f %f %f %f %f %f %f %f\n",
           p[indices[1]],
           width,
           x - d0mass,
           pUsingX,
           pUsingRealMass,
           bwFactor(pUsingRealMass),
           phspfactor,
           phspMassSq,
           phspGammaSq,
           ret);
    */
#endif

    //  if (gpuDebug & 1) printf("[%i, %i] KinLimitBW: %f %f %f %f %f\n", BLOCKIDX, THREADIDX, x, mean, width, d0mass, pimass, ret);
    return ret;
}

MEM_DEVICE device_function_ptr ptr_to_KinLimitBW = device_KinLimitBW;

__host__ KinLimitBWPdf::KinLimitBWPdf(std::string n, Variable* _x, Variable* mean, Variable* width)
    : GooPdf(_x, n) {
    registerParameter(mean);
    registerParameter(width);

    std::vector<unsigned int> pindices;
    pindices.push_back(mean->getIndex());
    pindices.push_back(width->getIndex());
    pindices.push_back(registerConstants(2));
    setMasses(1.8645, 0.13957);
    GET_FUNCTION_ADDR(ptr_to_KinLimitBW);
    initialise(pindices);
}

__host__ void KinLimitBWPdf::setMasses(fptype bigM, fptype smallM) {
    fptype constants[2];
    constants[0] = bigM;
    constants[1] = smallM;
    MEMCPY_TO_SYMBOL(functorConstants, constants, 2*sizeof(fptype), cIndex*sizeof(fptype), cudaMemcpyHostToDevice);
}
