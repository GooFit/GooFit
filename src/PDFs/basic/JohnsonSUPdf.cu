#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/basic/JohnsonSUPdf.h>

namespace GooFit {

const fptype SQRT2PI = 2.506628;

__device__ auto device_JohnsonSU(fptype *evt, ParameterContainer &pc) -> fptype {
    int id     = pc.getObservable(0);
    fptype x   = RO_CACHE(evt[id]);
    fptype _Jm = pc.getParameter(0);
    fptype _Js = pc.getParameter(1);
    fptype _Jg = pc.getParameter(2);
    fptype _Jd = pc.getParameter(3);
    pc.incrementIndex(1, 4, 0, 1, 1);

    fptype px       = (x - _Jm) / _Js;
    fptype px2      = px * px;
    fptype sqrt_arg = sqrt(1 + px2);
    fptype inv_sinh = log(px + sqrt_arg);
    fptype gaus_arg = _Jg + _Jd * inv_sinh;

    return _Jd / (_Js * SQRT2PI * sqrt_arg) * exp(-0.5 * gaus_arg * gaus_arg);
}

__device__ device_function_ptr ptr_to_JohnsonSU = device_JohnsonSU;

__host__
JohnsonSUPdf::JohnsonSUPdf(std::string n, Observable _x, Variable mean, Variable sigma, Variable gamma, Variable delta)
    : GooPdf("JohnsonSUPdf", n, _x, mean, sigma, gamma, delta) {
    registerFunction("ptr_to_JohnsonSU", ptr_to_JohnsonSU);

    initialize();
}

__host__ auto JohnsonSUPdf::integrate(fptype lo, fptype hi) const -> fptype {
    return 1.0; // Analytic integral included in device function! (Correct for minus to plus inf.)
}
} // namespace GooFit
