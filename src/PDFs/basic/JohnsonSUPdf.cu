#include <goofit/PDFs/basic/JohnsonSUPdf.h>
#include <goofit/PDFs/ParameterContainer.h>

namespace GooFit {

const fptype SQRT2PI = 2.506628;

__device__ fptype device_JohnsonSU(fptype *evt, ParameterContainer &pc) {
    int id     = pc.getObservable(0);
    fptype _Jm = pc.getParameter(0);
    fptype _Js = pc.getParameter(1);
    fptype _Jg = pc.getParameter(2);
    fptype _Jd = pc.getParameter(3);

    // we are using index 0.  If we need a different idx, we need to pass that information along.
    fptype x = evt[id];

    pc.incrementIndex(1, 4, 0, 1, 1);

    fptype px       = (x - _Jm) / _Js;
    fptype px2      = px * px;
    fptype sqrt_arg = sqrt(1 + px2);
    fptype inv_sinh = log(px + sqrt_arg);
    fptype gaus_arg = _Jg + _Jd * inv_sinh;
    // if ((gpuDebug & 1) && (0 == BLOCKIDX) && (0 == THREADIDX))
    // if (gpuDebug & 1)
    // printf("Johnson SU: %f %f %f %f | %f %f %i\n", _Jm, _Js, _Jg, _Jd, x, _Jd / (_Js * SQRT2PI * sqrt_arg) * exp(-0.5
    // * gaus_arg * gaus_arg), indices[2 + indices[0]]);
    // printf("Johnson SU: %f %f %f %f | %f %f %f %f\n", _Jm, _Js, _Jg, _Jd, x, _Jd / (_Js * SQRT2PI * sqrt_arg) *
    // exp(-0.5 * gaus_arg * gaus_arg), cudaArray[indices[1]], cudaArray[indices[2]]);
    return _Jd / (_Js * SQRT2PI * sqrt_arg) * exp(-0.5 * gaus_arg * gaus_arg);
}

__device__ device_function_ptr ptr_to_JohnsonSU = device_JohnsonSU;

__host__
JohnsonSUPdf::JohnsonSUPdf(std::string n, Observable _x, Variable mean, Variable sigma, Variable gamma, Variable delta)
    : GooPdf(n, _x) {
    registerParameter(mean);
    registerParameter(sigma);
    registerParameter(gamma);
    registerParameter(delta);

    initialize();
}

__host__ void JohnsonSUPdf::recursiveSetIndices() {
    GET_FUNCTION_ADDR(ptr_to_JohnsonSU);

    GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName(), "ptr_to_JohnsonSU");
    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx                               = num_device_functions++;

    populateArrays();
}

__host__ fptype JohnsonSUPdf::integrate(fptype lo, fptype hi) const {
    return 1.0; // Analytic integral included in device function! (Correct for minus to plus inf.)
}
} // namespace GooFit
