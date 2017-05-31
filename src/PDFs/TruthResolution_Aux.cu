#include "goofit/PDFs/physics/TruthResolution_Aux.h"

namespace GooFit {

__device__ fptype device_truth_resolution(fptype coshterm,
                                          fptype costerm,
                                          fptype sinhterm,
                                          fptype sinterm,
                                          fptype tau,
                                          fptype dtime,
                                          fptype xmixing,
                                          fptype ymixing,
                                          fptype /*sigma*/,
                                          fptype * /*p*/,
                                          unsigned int * /*indices*/) {
    fptype ret = 0;
    dtime /= tau;
    ret += coshterm * cosh(ymixing * dtime);
    ret += costerm * cos(xmixing * dtime);
    ret -= 2 * sinhterm * sinh(ymixing * dtime);
    ret -= 2 * sinterm
           * sin(xmixing * dtime); // Notice sign difference wrt to Mikhail's code, because I have AB* and he has A*B.
    ret *= exp(-dtime);

    // printf("device_truth_resolution %f %f %f %f %f\n", coshterm, costerm, sinhterm, sinterm, dtime);
    return ret;
}

__device__ fptype device_truth_resolution_average_tau(
    fptype A2, fptype B2, fptype ABr, fptype ABi, fptype xmixing, fptype ymixing, fptype tau) {
    fptype a          = A2 - B2;
    fptype b          = 2 * ABi;
    fptype c          = A2 + B2;
    fptype d          = 2 * ABr;
    fptype averagetau = ((xmixing * xmixing + 1) * (ymixing * ymixing - 1)
                         * (((a * tau * (xmixing * xmixing - 1.) + 2. * b * tau * xmixing))
                                / ((xmixing * xmixing + 1.) * (xmixing * xmixing + 1.))
                            + (c * (-(tau * ymixing * ymixing) - tau) + d * (2. * tau) * ymixing)
                                  / ((ymixing * ymixing - 1.) * (ymixing * ymixing - 1.))))
                        / ((ymixing * ymixing - 1) * (b * xmixing - a) + (xmixing * xmixing + 1) * (c - d * ymixing));
    // printf("device avg tau: %.5g with A2: %.5g, B2: %.5g, ABr:%.5g, ABi:%.5g, x:%.5g, y:%.5g, tau:%.5g \n",
    // averagetau, A2, B2, ABr, ABi, xmixing, ymixing, tau);
    return averagetau;
}

__device__ device_resfunction_ptr ptr_to_truth     = device_truth_resolution;
__device__ device_calc_tau_fcn_ptr ptr_to_calc_tau = device_truth_resolution_average_tau;

TruthResolution::TruthResolution()
    : MixingTimeResolution() {
    GET_FUNCTION_ADDR(ptr_to_truth);
    initIndex();
    GET_FUNCTION_ADDR(ptr_to_calc_tau);
    setCalcTauIdx(GooPdf::findFunctionIdx(host_fcn_ptr));
}
TruthResolution::~TruthResolution() = default;

fptype TruthResolution::normalisation(
    fptype di1, fptype di2, fptype di3, fptype di4, fptype tau, fptype xmixing, fptype ymixing) const {
    fptype timeIntegralOne = tau / (1 - ymixing * ymixing);
    fptype timeIntegralTwo = tau / (1 + xmixing * xmixing);
    fptype timeIntegralThr = ymixing * timeIntegralOne;
    fptype timeIntegralFou = xmixing * timeIntegralTwo;

    fptype ret = timeIntegralOne * (di1 + di2);
    ret += timeIntegralTwo * (di1 - di2);
    ret -= 2 * timeIntegralThr * di3;
    ret -= 2 * timeIntegralFou * di4;

    return ret;
}

} // namespace GooFit
