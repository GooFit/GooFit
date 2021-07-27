#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/basic/ExpGausPdf.h>

namespace GooFit {

__device__ auto device_ExpGaus(fptype *evt, ParameterContainer &pc) -> fptype {
    int id = pc.getObservable(0);

    fptype x     = RO_CACHE(evt[id]);
    fptype mean  = pc.getParameter(0);
    fptype sigma = pc.getParameter(1);
    fptype alpha = pc.getParameter(2);

    fptype ret    = 0.5 * alpha;
    fptype exparg = ret * (2 * mean + alpha * sigma * sigma - 2 * x);
    fptype erfarg = (mean + alpha * sigma * sigma - x) / (sigma * 1.4142135623);

    ret *= exp(exparg);
    ret *= erfc(erfarg);

    pc.incrementIndex(1, 3, 0, 1, 1);

    return ret;
}

__device__ device_function_ptr ptr_to_ExpGaus = device_ExpGaus;

ExpGausPdf::ExpGausPdf(std::string n, Observable _x, Variable mean, Variable sigma, Variable tau)
    : GooPdf("ExpGausPdf", n, _x, mean, sigma, tau) {
    registerFunction("ptr_to_ExpGaus", ptr_to_ExpGaus);

    initialize();
}

} // namespace GooFit
