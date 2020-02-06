#include <cmath>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/ThreeGaussResolution.h>

namespace GooFit {

const fptype R1o6 = 1.0 / 6.0;
#define SQRTPIo2 (1.0 / M_2_SQRTPI)
#define SQRT1o2PI (sqrt(0.5 * M_1_PI))

__device__ void gaussian(fptype &_P1,
                         fptype &_P2,
                         fptype &_P3,
                         fptype &_P4,
                         fptype _tau,
                         fptype adjTime,
                         fptype xmixing,
                         fptype ymixing,
                         fptype adjSigma,
                         fptype selbias) {
    fptype _1oSqrtA  = adjSigma * M_SQRT2;
    fptype _1oSigma  = 1 / adjSigma;
    fptype _1o2SqrtA = 0.5 * _1oSqrtA;
    fptype _1oSigma2 = _1oSigma * _1oSigma;
    fptype _NormG    = SQRT1o2PI * _1oSigma;

    fptype _C   = 0.5 * adjTime * adjTime * _1oSigma2;
    fptype _Bgn = -adjTime * _1oSigma2;

    fptype _Gamma = 1 / _tau;
    fptype _B     = _Gamma + _Bgn;

    fptype _u0  = _1o2SqrtA * (_B + selbias);
    fptype _u02 = _u0 * _u0;
    fptype _F   = _1oSqrtA * exp(-_C + _u02);

    fptype _Ig0 = SQRTPIo2 * erfc(_u0);
    fptype _Ig1 = 0.5 * exp(-_u02);
    fptype _Ig2 = _Ig1 * _u0 + 0.5 * _Ig0;
    fptype _Ig3 = _Ig1 * (_u02 + 1);

    fptype _R   = xmixing * _Gamma * _1oSqrtA;
    fptype _R2  = _R * _R;
    fptype _It0 = _F * _Ig0;
    fptype _It1 = _F * _R * (_Ig1 - _u0 * _Ig0);
    fptype _It2 = _F * _R2 * (_Ig2 - _u0 * _Ig1 * 2 + _u02 * _Ig0);
    fptype _It3 = _F * _R2 * _R * (_Ig3 - _u0 * _Ig2 * 3 + _u02 * _Ig1 * 3 - _u0 * _u02 * _Ig0);

    // fptype _P0   = _NormG *  _It0;
    _P2 = _NormG * (_It0 - 0.5 * _It2);
    _P4 = _NormG * (_It1 - R1o6 * _It3);

    fptype _u0py = _1o2SqrtA * (_B - ymixing * _Gamma + selbias);
    fptype _u0my = _1o2SqrtA * (_B + ymixing * _Gamma + selbias);
    fptype _Fpy  = _1oSqrtA * exp(-_C + _u0py * _u0py);
    fptype _Fmy  = _1oSqrtA * exp(-_C + _u0my * _u0my);
    fptype _Ipy  = _Fpy * SQRTPIo2 * erfc(_u0py);
    fptype _Imy  = _Fmy * SQRTPIo2 * erfc(_u0my);
    _P1          = _NormG * 0.5 * (_Ipy + _Imy);
    _P3          = _NormG * 0.5 * (_Ipy - _Imy);
}

__device__ fptype device_threegauss_resolution(fptype coshterm,
                                               fptype costerm,
                                               fptype sinhterm,
                                               fptype sinterm,
                                               fptype tau,
                                               fptype dtime,
                                               fptype xmixing,
                                               fptype ymixing,
                                               fptype sigma,
                                               ParameterContainer &pc) {
    fptype coreFraction    = pc.getParameter(0);
    fptype tailFraction    = (1 - coreFraction) * pc.getParameter(1);
    fptype outlFraction    = 1 - coreFraction - tailFraction;
    fptype coreBias        = pc.getParameter(2);
    fptype coreScaleFactor = pc.getParameter(3);
    fptype tailBias        = pc.getParameter(4);
    fptype tailScaleFactor = pc.getParameter(5);
    fptype outlBias        = pc.getParameter(6);
    fptype outlScaleFactor = pc.getParameter(7);
    fptype selbias         = pc.getParameter(8);

    fptype cp1 = 0;
    fptype cp2 = 0;
    fptype cp3 = 0;
    fptype cp4 = 0;
    gaussian(cp1, cp2, cp3, cp4, tau, dtime - coreBias * sigma, xmixing, ymixing, coreScaleFactor * sigma, selbias);
    fptype tp1 = 0;
    fptype tp2 = 0;
    fptype tp3 = 0;
    fptype tp4 = 0;
    gaussian(tp1, tp2, tp3, tp4, tau, dtime - tailBias * sigma, xmixing, ymixing, tailScaleFactor * sigma, selbias);
    fptype op1 = 0;
    fptype op2 = 0;
    fptype op3 = 0;
    fptype op4 = 0;
    gaussian(op1, op2, op3, op4, tau, dtime - outlBias * sigma, xmixing, ymixing, outlScaleFactor * sigma, selbias);

    fptype _P1 = coreFraction * cp1 + tailFraction * tp1 + outlFraction * op1;
    fptype _P2 = coreFraction * cp2 + tailFraction * tp2 + outlFraction * op2;
    fptype _P3 = coreFraction * cp3 + tailFraction * tp3 + outlFraction * op3;
    fptype _P4 = coreFraction * cp4 + tailFraction * tp4 + outlFraction * op4;

    fptype ret = 0;
    ret += coshterm * _P1;
    ret += costerm * _P2;
    ret -= 2 * sinhterm * _P3;
    ret -= 2 * sinterm * _P4; // Notice sign difference wrt to Mikhail's code, because I have AB* and he has A*B.

    // pc.incrementIndex (1, 8, 0, 0, 1);

    // cuPrintf("device_threegauss_resolution %f %f %f %f %f\n", coshterm, costerm, sinhterm, sinterm, dtime);
    return ret;
}

__device__ device_resfunction_ptr ptr_to_threegauss = device_threegauss_resolution;

ThreeGaussResolution::ThreeGaussResolution(
    Variable cf, Variable tf, Variable cb, Variable cs, Variable tb, Variable ts, Variable ob, Variable os, Variable sb)
    : MixingTimeResolution("ThreeGaussResolution", cf, tf, cb, cs, tb, ts, ob, os, sb)
    , selectionBias(sb) {
    initIndex();

    registerFunction("ptr_to_threegauss", ptr_to_threegauss);
}
ThreeGaussResolution::~ThreeGaussResolution() = default;

fptype ThreeGaussResolution::normalization(
    fptype di1, fptype di2, fptype di3, fptype di4, fptype tau, fptype xmixing, fptype ymixing) const {
    // NB! In thesis notation, A_1 = (A + B), A_2 = (A - B).
    // Here di1 = |A^2|, di2 = |B^2|, di3,4 = Re,Im(AB^*).
    // Distinction between numerical subscribts and A,B is crucial
    // for comparing thesis math to this math!

    // fptype timeIntegralOne = tau / (1 - ymixing * ymixing);
    // fptype timeIntegralTwo = tau / (1 + xmixing * xmixing);
    // fptype timeIntegralThr = ymixing * timeIntegralOne;
    // fptype timeIntegralFou = xmixing * timeIntegralTwo;

    fptype selBias = selectionBias.getValue();
    fptype timeIntegralOne
        = (selBias + 1 / tau) / (selBias * selBias + 2 * selBias / tau + (1 - ymixing * ymixing) / (tau * tau));
    fptype timeIntegralTwo
        = (selBias + 1 / tau) / (selBias * selBias + 2 * selBias / tau + (1 + xmixing * xmixing) / (tau * tau));
    fptype timeIntegralThr
        = (ymixing / tau) / (selBias * selBias + 2 * selBias / tau + (1 - ymixing * ymixing) / (tau * tau));
    fptype timeIntegralFou
        = (xmixing / tau) / (selBias * selBias + 2 * selBias / tau + (1 + xmixing * xmixing) / (tau * tau));

    fptype ret = timeIntegralOne * (di1 + di2); // ~ |A|^2 + |B|^2
    ret += timeIntegralTwo * (di1 - di2);       // ~ Re(A_1 A_2^*)
    ret -= 2 * timeIntegralThr * di3;           // ~ |A|^2 - |B|^2
    ret -= 2 * timeIntegralFou * di4;           // ~ Im(A_1 A_2^*)

    return ret;
}
} // namespace GooFit
