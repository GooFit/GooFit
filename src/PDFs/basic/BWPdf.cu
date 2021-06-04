#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/basic/BWPdf.h>

namespace GooFit {

__device__ auto device_BW(fptype *evt, ParameterContainer &pc) -> fptype {
    int id = pc.getObservable(0);

    fptype x      = RO_CACHE(evt[id]);
    fptype mean   = pc.getParameter(0);
    fptype gamma  = pc.getParameter(1);
    fptype rootPi = -2. * atan2(-1.0, 0.0);
    fptype ret    = (gamma / ((x - mean) * (x - mean) + gamma * gamma / 4)) / (2 * rootPi);

    pc.incrementIndex(1, 2, 0, 1, 1);

    return ret;
}

__device__ device_function_ptr ptr_to_BW = device_BW;

__host__ BWPdf::BWPdf(std::string n, Observable _x, Variable mean, Variable width)
    : GooPdf("BWPdf", n, _x, mean, width) {
    registerFunction("ptr_to_BW", ptr_to_BW);

    initialize();
}

} // namespace GooFit
