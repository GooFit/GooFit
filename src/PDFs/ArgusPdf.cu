#include "goofit/PDFs/basic/ArgusPdf.h"
#include "goofit/Variable.h"

namespace GooFit {

__device__ fptype device_Argus_Upper(fptype *evt, fptype *p, unsigned int *indices) {
    fptype x  = evt[indices[2 + indices[0]]];
    fptype m0 = p[indices[1]];

    double t = x / m0;

    if(t >= 1)
        return 0;

    fptype slope = p[indices[2]];
    fptype power = p[indices[3]];
    t            = 1 - t * t;
    // printf("device_Argus_Upper %f %f %f %f %f\n", x, m0, slope, t, x * pow(t, power) * exp(slope * t));

    return x * pow(t, power) * exp(slope * t);
}

__device__ fptype device_Argus_Lower(fptype *evt, fptype *p, unsigned int *indices) {
    fptype x  = evt[indices[2 + indices[0]]];
    fptype m0 = p[indices[1]];

    // printf("Argus: %i %i %f %f\n", indices[0], indices[2 + indices[0]], x, m0);
    // printf("Argus: %i %i\n", indices[0], indices[2 + indices[0]]);
    // return 1;

    fptype t = x / m0;

    if(t <= 1)
        return 0;

    t *= t;
    t -= 1;

    fptype slope = p[indices[2]];
    fptype power = p[indices[3]];
    fptype ret   = x * pow(t, power) * exp(slope * t);
    // if ((0 == THREADIDX) && (0 == BLOCKIDX) && (callnumber < 1)) cuPrintf("device_Argus_Lower %i %i %f %f %f %f
    // %f\n", indices[1], indices[2], x, m0, slope, t, ret);
    // if (isnan(ret)) printf("NaN Argus: %f %f %f %f %f %f %f\n", x, m0, t, slope, power, pow(t, power), exp(slope*t));
    // if ((0 == THREADIDX) && (0 == BLOCKIDX) && (gpuDebug & 1))
    // printf("(%i, %i) device_Argus_Lower %f %f %f %f %f\n", BLOCKIDX, THREADIDX, x, m0, slope, t, x * pow(t, power) *
    // exp(slope * t));

    return ret;
}

__device__ device_function_ptr ptr_to_Argus_Upper = device_Argus_Upper;
__device__ device_function_ptr ptr_to_Argus_Lower = device_Argus_Lower;

__host__ ArgusPdf::ArgusPdf(std::string n, Variable *_x, Variable *m0, Variable *slope, bool upper, Variable *power)
    : GooPdf(_x, n) {
    registerParameter(m0);
    registerParameter(slope);

    if(!power)
        power = new Variable(n + "powervar", 0.5);

    registerParameter(power);

    std::vector<unsigned int> pindices;
    pindices.push_back(m0->getIndex());
    pindices.push_back(slope->getIndex());
    pindices.push_back(power->getIndex());

    if(upper)
        GET_FUNCTION_ADDR(ptr_to_Argus_Upper);
    else
        GET_FUNCTION_ADDR(ptr_to_Argus_Lower);

    initialize(pindices);
}

fptype argus_lower_helper(fptype x, fptype m0, fptype slope, fptype power) {
    fptype t = x / m0;

    if(t <= 1)
        return 0;

    t *= t;
    t -= 1;

    fptype ret = x * pow(t, power) * exp(slope * t);

    return ret;
}

__host__ double ArgusPdf::integrate(fptype lo, fptype hi) const {
    double norm           = 0;
    unsigned int *indices = host_indices + parameters;
    fptype m0             = host_params[indices[1]];
    fptype slope          = host_params[indices[2]];
    fptype power          = host_params[indices[3]];

    for(int j = 0; j < integrationBins; ++j) {
        double x = hi;
        x -= lo;
        x /= integrationBins;
        x *= j;
        x += lo;
        norm += argus_lower_helper(x, m0, slope, power);
    }

    norm *= ((hi - lo) / integrationBins);
    return norm;
}
} // namespace GooFit
