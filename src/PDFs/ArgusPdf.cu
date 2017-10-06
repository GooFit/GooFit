#include "goofit/PDFs/basic/ArgusPdf.h"
#include "goofit/Variable.h"

namespace GooFit {


__device__ fptype device_Argus_Upper(fptype* evt, ParameterContainer &pc) {
    int id = pc.observables[pc.observableIdx + 1];

    fptype x = evt[id];
    fptype m0 = pc.parameters[pc.parameterIdx + 1];

    double t = x / m0;

    fptype slope = pc.parameters[pc.parameterIdx + 2];
    fptype power = pc.parameters[pc.parameterIdx + 3];

    pc.incrementIndex (1, 3, 0, 1, 1);

    if(t >= 1)
        return 0;

    t = 1 - t*t;

    return x * pow(t, power) * exp(slope * t);
}

__device__ fptype device_Argus_Lower(fptype* evt, ParameterContainer &pc) {
    int id = pc.observables[pc.observableIdx + 1];

    fptype x = evt[id];
    fptype m0 = pc.parameters[pc.parameterIdx + 1];

    fptype t = x / m0;

    fptype slope = pc.parameters[pc.parameterIdx + 2];
    fptype power = pc.parameters[pc.parameterIdx + 3];

    pc.incrementIndex (1, 3, 0, 1, 1);

    if(t <= 1)
        return 0;

    t *= t;
    t -= 1;

    fptype ret = x * pow(t, power) * exp(slope * t);

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

    if(upper) {
        ArgusType = 1;
        GET_FUNCTION_ADDR(ptr_to_Argus_Upper);
    }
    else {
        ArgusType = 0;
        GET_FUNCTION_ADDR(ptr_to_Argus_Lower);
    }

    initialize(pindices);
}

__host__ void ArgusPdf::recursiveSetIndices () {
    if (ArgusType == 1) {
        GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName (), "ptr_to_Argus_Upper");
        GET_FUNCTION_ADDR(ptr_to_Argus_Upper);
    } else if (ArgusType == 0) {
        GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName (), "ptr_to_Argus_Lower");
        GET_FUNCTION_ADDR(ptr_to_Argus_Lower);
    }

    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx = num_device_functions ++;

    populateArrays ();
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
    double norm = 0;
    fptype m0    = host_parameters[parametersIdx + 1];
    fptype slope = host_parameters[parametersIdx + 2];
    fptype power = host_parameters[parametersIdx + 3];

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
