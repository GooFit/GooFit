#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/basic/ArgusPdf.h>
#include <goofit/Variable.h>

namespace GooFit {

__device__ auto device_Argus_Upper(fptype *evt, ParameterContainer &pc) -> fptype {
    int id = pc.getObservable(0);

    fptype x     = RO_CACHE(evt[id]);
    fptype m0    = pc.getParameter(0);
    fptype slope = pc.getParameter(1);
    fptype power = pc.getParameter(2);

    pc.incrementIndex(1, 3, 0, 1, 1);

    double t = x / m0;

    if(t >= 1)
        return 0;

    t = 1 - t * t;

    return x * pow(t, power) * exp(slope * t);
}

__device__ auto device_Argus_Lower(fptype *evt, ParameterContainer &pc) -> fptype {
    int id = pc.getObservable(0);

    fptype x     = RO_CACHE(evt[id]);
    fptype m0    = pc.getParameter(0);
    fptype slope = pc.getParameter(1);
    fptype power = pc.getParameter(2);

    pc.incrementIndex(1, 3, 0, 1, 1);

    fptype t = x / m0;
    if(t <= 1)
        return 0;

    t *= t;
    t -= 1;

    fptype ret = x * pow(t, power) * exp(slope * t);

    return ret;
}

__device__ device_function_ptr ptr_to_Argus_Upper = device_Argus_Upper;
__device__ device_function_ptr ptr_to_Argus_Lower = device_Argus_Lower;

__host__ ArgusPdf::ArgusPdf(std::string n, Observable _x, Variable m0, Variable slope, bool upper)
    : ArgusPdf(n, _x, m0, slope, upper, Variable(n + "powervar", 0.5)) {}

__host__ ArgusPdf::ArgusPdf(std::string n, Observable _x, Variable m0, Variable slope, bool upper, Variable power)
    : GooPdf("ArgusPdf", n, _x, m0, slope, power) {
    if(upper)
        registerFunction("ptr_to_Argus_Upper", ptr_to_Argus_Upper);
    else
        registerFunction("ptr_to_Argus_Lower", ptr_to_Argus_Lower);

    initialize();
}

auto argus_lower_helper(fptype x, fptype m0, fptype slope, fptype power) -> fptype {
    fptype t = x / m0;

    if(t <= 1)
        return 0;

    t *= t;
    t -= 1;

    fptype ret = x * pow(t, power) * exp(slope * t);

    return ret;
}

__host__ auto ArgusPdf::integrate(fptype lo, fptype hi) const -> double {
    double norm  = 0;
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
