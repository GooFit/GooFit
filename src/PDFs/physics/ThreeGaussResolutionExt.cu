#include <cmath>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/ThreeGaussResolutionExt.h>

namespace GooFit {

const fptype R1o6 = 1.0 / 6.0;
#define SQRTPIo2 (1.0 / M_2_SQRTPI)
#define SQRT1o2PI (sqrt(0.5 * M_1_PI))

__device__ void gaussian2(fptype &_P1,
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

// convolution with constant decay acceptance in region 0<t<T_thres
__device__ void gaussian_low(fptype &_P1,
                             fptype &_P2,
                             fptype &_P3,
                             fptype &_P4,
                             fptype _tau,
                             fptype adjTime,
                             fptype xmixing,
                             fptype ymixing,
                             fptype adjSigma,
                             fptype selbias,
                             fptype Tthres,
                             fptype C) {
    fptype preConst  = C * exp(selbias * Tthres);
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
    fptype u1   = 1. / _1oSqrtA * Tthres + _u0;
    fptype _u02 = _u0 * _u0;
    fptype _F   = _1oSqrtA * exp(-_C + _u02);

    fptype _Ig0 = SQRTPIo2 * erfc(_u0) - SQRTPIo2 * erfc(u1);
    fptype _Ig1 = 0.5 * exp(-_u02) - 0.5 * exp(-u1 * u1);
    fptype _Ig2
        = 0.5 * (_u0 * exp(-1 * _u02) + SQRTPIo2 * erfc(_u0)) - 0.5 * (u1 * exp(-u1 * u1) + SQRTPIo2 * erfc(u1));
    fptype _Ig3 = 0.5 * exp(-_u02) * (_u02 + 1) - 0.5 * exp(-u1 * u1) * (u1 * u1 + 1);

    fptype _R   = xmixing * _Gamma * _1oSqrtA;
    fptype _R2  = _R * _R;
    fptype _It0 = _F * _Ig0;
    fptype _It1 = _F * _R * (_Ig1 - _u0 * _Ig0);
    fptype _It2 = _F * _R2 * (_Ig2 - _u0 * _Ig1 * 2 + _u02 * _Ig0);
    fptype _It3 = _F * _R2 * _R * (_Ig3 - _u0 * _Ig2 * 3 + _u02 * _Ig1 * 3 - _u0 * _u02 * _Ig0);

    // fptype _P0   = _NormG *  _It0;
    _P2 = preConst * _NormG * (_It0 - 0.5 * _It2);
    _P4 = preConst * _NormG * (_It1 - R1o6 * _It3);

    fptype _u0py = _1o2SqrtA * (_B - ymixing * _Gamma + selbias);
    fptype _u0my = _1o2SqrtA * (_B + ymixing * _Gamma + selbias);
    fptype u1py  = _u0py + 1. / _1oSqrtA * Tthres;
    fptype u1my  = _u0my + 1. / _1oSqrtA * Tthres;
    fptype _Fpy  = _1oSqrtA * exp(-_C + _u0py * _u0py);
    fptype _Fmy  = _1oSqrtA * exp(-_C + _u0my * _u0my);
    fptype _Ipy  = _Fpy * SQRTPIo2 * (erfc(_u0py) - erfc(u1py));
    fptype _Imy  = _Fmy * SQRTPIo2 * (erfc(_u0my) - erfc(u1my));
    _P1          = preConst * _NormG * 0.5 * (_Ipy + _Imy);
    _P3          = preConst * _NormG * 0.5 * (_Ipy - _Imy);
}

// convolution with exponential decay acceptance in region T_thres<t<inf
__device__ void gaussian_high(fptype &_P1,
                              fptype &_P2,
                              fptype &_P3,
                              fptype &_P4,
                              fptype _tau,
                              fptype adjTime,
                              fptype xmixing,
                              fptype ymixing,
                              fptype adjSigma,
                              fptype selbias,
                              fptype Tthres,
                              fptype C) {
    fptype preConst  = C * exp(selbias * Tthres);
    fptype _1oSqrtA  = adjSigma * M_SQRT2; // 1/sqrt(r)
    fptype _1oSigma  = 1 / adjSigma;
    fptype _1o2SqrtA = 0.5 * _1oSqrtA; // 1/(2 sqrt(r))
    fptype _1oSigma2 = _1oSigma * _1oSigma;
    fptype _NormG    = SQRT1o2PI * _1oSigma;

    fptype _C   = 0.5 * adjTime * adjTime * _1oSigma2;
    fptype _Bgn = -adjTime * _1oSigma2;

    fptype _Gamma = 1 / _tau;
    fptype _B     = _Gamma + _Bgn;

    fptype _u0  = _1o2SqrtA * (_B + selbias);
    fptype u1   = 1. / _1oSqrtA * Tthres + _u0;
    fptype _u02 = _u0 * _u0;
    fptype _F   = _1oSqrtA * exp(-_C + _u02);

    fptype _Ig0 = SQRTPIo2 * erfc(u1);
    fptype _Ig1 = 0.5 * exp(-1 * u1 * u1);
    fptype _Ig2 = _Ig1 * 1 + 0.5 * _Ig0;
    fptype _Ig3 = _Ig1 * (u1 + 1);

    fptype _R   = xmixing * _Gamma * _1oSqrtA;
    fptype _R2  = _R * _R;
    fptype _It0 = _F * _Ig0;
    fptype _It1 = _F * _R * (_Ig1 - _u0 * _Ig0);
    fptype _It2 = _F * _R2 * (_Ig2 - _u0 * _Ig1 * 2 + _u02 * _Ig0);
    fptype _It3 = _F * _R2 * _R * (_Ig3 - _u0 * _Ig2 * 3 + _u02 * _Ig1 * 3 - _u0 * _u02 * _Ig0);

    // fptype _P0   = _NormG *  _It0;
    _P2 = preConst * _NormG * (_It0 - 0.5 * _It2);
    _P4 = preConst * _NormG * (_It1 - R1o6 * _It3);

    fptype _u0py = _1o2SqrtA * (_B - ymixing * _Gamma + selbias);
    fptype _u0my = _1o2SqrtA * (_B + ymixing * _Gamma + selbias);
    fptype u1py  = _u0py + 1. / _1oSqrtA * Tthres;
    fptype u1my  = _u0my + 1. / _1oSqrtA * Tthres;
    fptype _Fpy  = _1oSqrtA * exp(-_C + _u0py * _u0py);
    fptype _Fmy  = _1oSqrtA * exp(-_C + _u0my * _u0my);
    fptype _Ipy  = _Fpy * SQRTPIo2 * erfc(u1py);
    fptype _Imy  = _Fmy * SQRTPIo2 * erfc(u1my);
    _P1          = preConst * _NormG * 0.5 * (_Ipy + _Imy);
    _P3          = preConst * _NormG * 0.5 * (_Ipy - _Imy);
}

__device__ fptype device_threegauss_resolutionext(fptype coshterm,
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
    fptype tailFraction    = pc.getParameter(1);
    fptype outlFraction    = 1 - coreFraction - tailFraction;
    fptype coreBias        = pc.getParameter(2);
    fptype coreScaleFactor = pc.getParameter(3);
    fptype tailBias        = pc.getParameter(4);
    fptype tailScaleFactor = pc.getParameter(5);
    fptype outlBias        = pc.getParameter(6);
    fptype outlScaleFactor = pc.getParameter(7);
    fptype selbias_low     = pc.getParameter(8);
    fptype selbias_high    = pc.getParameter(9);
    fptype Tthreshold      = pc.getParameter(10);
    fptype constantC       = pc.getParameter(11);

    fptype cp1_low = 0;
    fptype cp2_low = 0;
    fptype cp3_low = 0;
    fptype cp4_low = 0;
    fptype tp1_low = 0;
    fptype tp2_low = 0;
    fptype tp3_low = 0;
    fptype tp4_low = 0;
    fptype op1_low = 0;
    fptype op2_low = 0;
    fptype op3_low = 0;
    fptype op4_low = 0;

    fptype cp1_high = 0;
    fptype cp2_high = 0;
    fptype cp3_high = 0;
    fptype cp4_high = 0;
    fptype tp1_high = 0;
    fptype tp2_high = 0;
    fptype tp3_high = 0;
    fptype tp4_high = 0;
    fptype op1_high = 0;
    fptype op2_high = 0;
    fptype op3_high = 0;
    fptype op4_high = 0;

    gaussian_low(cp1_low,
                 cp2_low,
                 cp3_low,
                 cp4_low,
                 tau,
                 dtime - coreBias * sigma,
                 xmixing,
                 ymixing,
                 coreScaleFactor * sigma,
                 selbias_low,
                 Tthreshold,
                 constantC);
    gaussian_low(tp1_low,
                 tp2_low,
                 tp3_low,
                 tp4_low,
                 tau,
                 dtime - tailBias * sigma,
                 xmixing,
                 ymixing,
                 tailScaleFactor * sigma,
                 selbias_low,
                 Tthreshold,
                 constantC);
    gaussian_low(op1_low,
                 op2_low,
                 op3_low,
                 op4_low,
                 tau,
                 dtime - outlBias * sigma,
                 xmixing,
                 ymixing,
                 outlScaleFactor * sigma,
                 selbias_low,
                 Tthreshold,
                 constantC);

    gaussian_high(cp1_high,
                  cp2_high,
                  cp3_high,
                  cp4_high,
                  tau,
                  dtime - coreBias * sigma,
                  xmixing,
                  ymixing,
                  coreScaleFactor * sigma,
                  selbias_high,
                  Tthreshold,
                  constantC);
    gaussian_high(tp1_high,
                  tp2_high,
                  tp3_high,
                  tp4_high,
                  tau,
                  dtime - tailBias * sigma,
                  xmixing,
                  ymixing,
                  tailScaleFactor * sigma,
                  selbias_high,
                  Tthreshold,
                  constantC);
    gaussian_high(op1_high,
                  op2_high,
                  op3_high,
                  op4_high,
                  tau,
                  dtime - outlBias * sigma,
                  xmixing,
                  ymixing,
                  outlScaleFactor * sigma,
                  selbias_high,
                  Tthreshold,
                  constantC);

    fptype _P1 = coreFraction * (cp1_low + cp1_high) + tailFraction * (tp1_low + tp1_high)
                 + outlFraction * (op1_low + op1_high);
    fptype _P2 = coreFraction * (cp2_low + cp2_high) + tailFraction * (tp2_low + tp2_high)
                 + outlFraction * (op2_low + op2_high);
    fptype _P3 = coreFraction * (cp3_low + cp3_high) + tailFraction * (tp3_low + tp3_high)
                 + outlFraction * (op3_low + op3_high);
    fptype _P4 = coreFraction * (cp4_low + cp4_high) + tailFraction * (tp4_low + tp4_high)
                 + outlFraction * (op4_low + op4_high);

    fptype ret = 0;
    ret += coshterm * _P1;
    ret += costerm * _P2;
    ret -= 2 * sinhterm * _P3;
    ret -= 2 * sinterm * _P4; // Notice sign difference wrt to Mikhail's code, because I have AB* and he has A*B.

    // pc.incrementIndex (1, 8, 0, 0, 1);

    // cuPrintf("device_threegauss_resolution %f %f %f %f %f\n", coshterm, costerm, sinhterm, sinterm, dtime);
    return ret;
}

__device__ device_resfunction_ptr ptr_to_threegaussext = device_threegauss_resolutionext;

ThreeGaussResolutionExt::ThreeGaussResolutionExt(Variable cf,
                                                 Variable tf,
                                                 Variable cb,
                                                 Variable cs,
                                                 Variable tb,
                                                 Variable ts,
                                                 Variable ob,
                                                 Variable os,
                                                 Variable sb_low,
                                                 Variable sb_high,
                                                 Variable Tthres,
                                                 Variable constantC)
    : MixingTimeResolution(
        "ThreeGaussResolutionExt", cf, tf, cb, cs, tb, ts, ob, os, sb_low, sb_high, Tthres, constantC)
    , selectionBias_low(sb_low)
    , selectionBias_high(sb_high)
    , mTthreshold(Tthres)
    , mConstantC(constantC) {
    initIndex();

    registerFunction("ptr_to_threegaussext", ptr_to_threegaussext);
}
ThreeGaussResolutionExt::~ThreeGaussResolutionExt() = default;

fptype ThreeGaussResolutionExt::normalization(
    fptype di1, fptype di2, fptype di3, fptype di4, fptype tau, fptype xmixing, fptype ymixing) const {
    // NB! In thesis notation, A_1 = (A + B), A_2 = (A - B).
    // Here di1 = |A^2|, di2 = |B^2|, di3,4 = Re,Im(AB^*).
    // Distinction between numerical subscribts and A,B is crucial
    // for comparing thesis math to this math!

    // fptype timeIntegralOne = tau / (1 - ymixing * ymixing);
    // fptype timeIntegralTwo = tau / (1 + xmixing * xmixing);
    // fptype timeIntegralThr = ymixing * timeIntegralOne;
    // fptype timeIntegralFou = xmixing * timeIntegralTwo;

    fptype selBias_low   = selectionBias_low.getValue();
    fptype selBias_high  = selectionBias_high.getValue();
    fptype Tthres        = mTthreshold.getValue();
    fptype C             = mConstantC.getValue();
    fptype preConst_low  = C * exp(selBias_low * Tthres);
    fptype preConst_high = C * exp(selBias_high * Tthres);

    fptype Gamma         = 1. / tau;
    fptype gammaPlusBias = Gamma + selBias_low;

    fptype timeIntegralOne_low
        = 0.5 * preConst_low
          * (1. / (ymixing * Gamma - Gamma - selBias_low) * (exp(Tthres * (ymixing * Gamma - Gamma - selBias_low)) - 1)
             + 1. / (-ymixing * Gamma - Gamma - selBias_low)
                   * (exp(Tthres * (-ymixing * Gamma - Gamma - selBias_low)) - 1));

    fptype timeIntegralThr_low
        = 0.5 * preConst_low
          * (1. / (ymixing * Gamma - Gamma - selBias_low) * (exp(Tthres * (ymixing * Gamma - Gamma - selBias_low)) - 1)
             - 1. / (-ymixing * Gamma - Gamma - selBias_low)
                   * (exp(Tthres * (-ymixing * Gamma - Gamma - selBias_low)) - 1));

    fptype timeIntegralTwo_low
        = preConst_low
          * (exp(-gammaPlusBias * Tthres) / (gammaPlusBias * gammaPlusBias + xmixing * xmixing * Gamma * Gamma)
                 * (-gammaPlusBias * cos(xmixing * Gamma * Tthres) + xmixing * Gamma * sin(xmixing * Gamma * Tthres))
             - (-gammaPlusBias) / (gammaPlusBias * gammaPlusBias + xmixing * xmixing * Gamma * Gamma));

    fptype timeIntegralFour_low
        = preConst_low
          * (exp(-gammaPlusBias * Tthres) / (gammaPlusBias * gammaPlusBias + xmixing * xmixing * Gamma * Gamma)
                 * (-gammaPlusBias * sin(xmixing * Gamma * Tthres) - xmixing * Gamma * cos(xmixing * Gamma * Tthres))
             - (-xmixing * Gamma) / (gammaPlusBias * gammaPlusBias + xmixing * xmixing * Gamma * Gamma));

    gammaPlusBias = Gamma + selBias_high;

    fptype timeIntegralOne_high
        = 0.5 * preConst_high
          * (-1. / (ymixing * Gamma - Gamma - selBias_high) * (exp(Tthres * (ymixing * Gamma - Gamma - selBias_high)))
             - 1. / (-ymixing * Gamma - Gamma - selBias_high)
                   * (exp(Tthres * (-ymixing * Gamma - Gamma - selBias_high))));

    fptype timeIntegralThr_high
        = 0.5 * preConst_high
          * (-1. / (ymixing * Gamma - Gamma - selBias_high) * (exp(Tthres * (ymixing * Gamma - Gamma - selBias_high)))
             + 1. / (-ymixing * Gamma - Gamma - selBias_high)
                   * (exp(Tthres * (-ymixing * Gamma - Gamma - selBias_high))));

    fptype timeIntegralTwo_high
        = preConst_high
          * (-exp(-gammaPlusBias * Tthres) / (gammaPlusBias * gammaPlusBias + xmixing * xmixing * Gamma * Gamma)
             * (-gammaPlusBias * cos(xmixing * Gamma * Tthres) + xmixing * Gamma * sin(xmixing * Gamma * Tthres)));

    fptype timeIntegralFour_high
        = preConst_high
          * (-exp(-gammaPlusBias * Tthres) / (gammaPlusBias * gammaPlusBias + xmixing * xmixing * Gamma * Gamma)
             * (-gammaPlusBias * sin(xmixing * Gamma * Tthres) - xmixing * Gamma * cos(xmixing * Gamma * Tthres)));

    fptype timeIntegralOne = timeIntegralOne_low + timeIntegralOne_high;
    fptype timeIntegralTwo = timeIntegralTwo_low + timeIntegralTwo_high;
    fptype timeIntegralThr = timeIntegralThr_low + timeIntegralThr_high;
    fptype timeIntegralFou = timeIntegralFour_low + timeIntegralFour_high;

    fptype ret = timeIntegralOne * (di1 + di2); // ~ |A|^2 + |B|^2
    ret += timeIntegralTwo * (di1 - di2);       // ~ |A|^2 - |B|^2
    ret -= 2 * timeIntegralThr * di3;           // ~ Re(A_1 A_2^*)
    ret -= 2 * timeIntegralFou * di4;           // ~ Im(A_1 A_2^*)

    return ret;
}
} // namespace GooFit
