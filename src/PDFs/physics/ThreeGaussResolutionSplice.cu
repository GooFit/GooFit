#include <cmath>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/ThreeGaussResolutionSplice.h>

namespace GooFit {

const fptype R1o6 = 1.0 / 6.0;
#define SQRTPIo2 (1.0 / M_2_SQRTPI)
#define SQRT1o2PI (sqrt(0.5 * M_1_PI))
#define SQRT1oPI (sqrt(M_1_PI))

/*
TODO: make sure that acceptance function goes to zero
*/

__device__ fptype BigPhi(fptype u0, fptype u1) {
    fptype res = 0.5 * (erfc(u0 / sqrt(2.)) - erfc(u1 / sqrt(2.)));
    return res;
}

__device__ fptype BigPhi(fptype x) {
    fptype res = 0.5 * (1. + erf(x / sqrt(2.)));
    return res;
}

__device__ fptype SmallPhi(fptype u0, fptype u1) {
    fptype res = SQRT1o2PI * (exp(-0.5 * u1 * u1) - exp(-0.5 * u0 * u0));
    return res;
}

__device__ fptype SmallPhi(fptype x) {
    fptype res = SQRT1o2PI * (exp(-0.5 * x * x));
    return res;
}

__device__ fptype DoubleFactorial(int j) {
    // printf("DoubleFactorial start \n");
    fptype res = 1.;
    while(j > 0) {
        res = res * j;
        j   = j - 2;
    }
    // printf("DoubleFactorial finish \n");
    return (fptype)res;
}

__device__ fptype FactorialOdd(int n, fptype u0, fptype u1) {
    // evaluate odd moments of Gaussian
    // printf("FactorialOdd start \n");
    int k                 = (n - 1) / 2;
    fptype constFactorial = DoubleFactorial(2 * k);
    fptype res            = 0.;
    for(int j = 0; j <= k; j++) {
        res = res + 1. / DoubleFactorial(2 * j) * (-SmallPhi(u1) * pow(u1, 2 * j) + SmallPhi(u0) * pow(u0, 2 * j));
    }
    res = res * constFactorial;
    // printf("FactorialOdd finish \n");
    return res;
}

__device__ fptype FactorialEven(int n, fptype u0, fptype u1) {
    // evaluate even moments of Gaussian
    // printf("FactorialOdd start \n");
    int k                 = (n - 2) / 2;
    fptype constFactorial = DoubleFactorial(2 * k + 1);
    fptype res            = 0.;
    for(int j = 0; j <= k; j++) {
        res = res
              + 1. / DoubleFactorial(2 * j + 1)
                    * (-SmallPhi(u1) * pow(u1, 2 * j + 1) + SmallPhi(u0) * pow(u0, 2 * j + 1));
    }

    res = res + BigPhi(u1) - BigPhi(u0);
    res = res * constFactorial;
    // printf("FactorialOdd finish \n");
    return res;
}

__device__ fptype EvaluateConvo(int n, fptype u0, fptype u1) {
    // evaluate moments of Gaussian function
    if(n == 0)
        return BigPhi(u1) - BigPhi(u0);
    if(n == 1)
        return -SmallPhi(u1) + SmallPhi(u0);
    if(n % 2 == 0)
        return FactorialEven(n, u0, u1);
    else
        return FactorialOdd(n, u0, u1);
}

__host__ __device__ fptype Factorial(int j) {
    if(j == 0)
        return 1.;
    // printf("evaluating factorial of %d \n", j);
    fptype res = 1.;
    while(j > 0) {
        res = res * j;
        // printf("factorial res is %d \n", res);
        j--;
    }
    return (fptype)res;
}

__device__ fptype BinomCoeff(int n, int k) {
    fptype res = Factorial(n);
    res        = res / (Factorial(k) * Factorial(n - k));
    return res;
}

__device__ fptype EvaluateAcceptance(int n, fptype _sqrt_r, fptype *_evaluatedConvo, fptype *_evaluatedPowU0) {
    fptype res = 0.;
    for(int k = 0; k <= n; k++) {
        res += BinomCoeff(n, k) * _evaluatedConvo[n - k] * _evaluatedPowU0[k];
    }
    res = res * pow(1. / _sqrt_r, n);
    return res;
}

__device__ void EvaluateKnot(fptype &_P1,
                             fptype &_P3,
                             fptype Gamma,
                             fptype knot_low,
                             fptype knot_high,
                             fptype a0,
                             fptype a1,
                             fptype a2,
                             fptype a3,
                             fptype tprime,
                             fptype sigma,
                             fptype y) {
    // sqrt_r = 1/sigma
    fptype sqrt_r = 1. / sigma;
    // r = 1/sigma^2
    fptype r        = sqrt_r * sqrt_r;
    fptype p_plus   = 2 * Gamma * (1. - y) - 2 * (tprime)*r;
    fptype p_minus  = 2 * Gamma * (1. + y) - 2 * (tprime)*r;
    fptype q        = (tprime) * (tprime)*r;
    fptype u0_plus  = p_plus / (2. * sqrt_r);
    fptype u0_minus = p_minus / (2. * sqrt_r);
    // factor 1/sqrt(2PI) absorbed in small phi?
    fptype preFactor       = 0.5 / sqrt_r / sigma;
    fptype preFactor_plus  = exp(-0.5 * (q - u0_plus * u0_plus));
    fptype preFactor_minus = exp(-0.5 * (q - u0_minus * u0_minus));

    fptype lowLim  = sqrt_r * knot_low + u0_plus;
    fptype highLim = sqrt_r * knot_high + u0_plus;
    fptype evaluatedConvo[4];
    fptype evaluatedPowU0[4];
    for(int i = 0; i < 4; i++) {
        evaluatedConvo[i] = EvaluateConvo(i, lowLim, highLim);
    }
    evaluatedPowU0[0] = 1.;
    for(int i = 1; i < 4; i++) {
        evaluatedPowU0[i] = evaluatedPowU0[i - 1] * -u0_plus;
    }

    /*
    fptype evaluatedAcceptance[4];
    evaluatedAcceptance[0] = BinomCoeff(0, 0) * EvaluateConvo(0, sqrt_r * knot_low+  u0_plus, _sqrt_r * knot_high +
    u0_plus) * pow(- u0_plus, 0); for(int i = 1; i<4; i++) { evaluatedAcceptance[i] = evaluatedAcceptance[i-1] +
    BinomCoeff(n, k) * EvaluateConvo(n - k, _sqrt_r * tlow + u0, _sqrt_r * thigh + u0) * pow(-u0, k);
    }
    */

    fptype Ipy_0 = EvaluateAcceptance(0, sqrt_r, evaluatedConvo, evaluatedPowU0);
    fptype Ipy_1 = EvaluateAcceptance(1, sqrt_r, evaluatedConvo, evaluatedPowU0);
    fptype Ipy_2 = EvaluateAcceptance(2, sqrt_r, evaluatedConvo, evaluatedPowU0);
    fptype Ipy_3 = EvaluateAcceptance(3, sqrt_r, evaluatedConvo, evaluatedPowU0);

    lowLim  = sqrt_r * knot_low + u0_minus;
    highLim = sqrt_r * knot_high + u0_minus;
    for(int i = 0; i < 4; i++) {
        evaluatedConvo[i] = EvaluateConvo(i, lowLim, highLim);
    }
    evaluatedPowU0[0] = 1.;
    for(int i = 1; i < 4; i++) {
        evaluatedPowU0[i] = evaluatedPowU0[i - 1] * -u0_minus;
    }

    fptype Imy_0 = EvaluateAcceptance(0, sqrt_r, evaluatedConvo, evaluatedPowU0);
    fptype Imy_1 = EvaluateAcceptance(1, sqrt_r, evaluatedConvo, evaluatedPowU0);
    fptype Imy_2 = EvaluateAcceptance(2, sqrt_r, evaluatedConvo, evaluatedPowU0);
    fptype Imy_3 = EvaluateAcceptance(3, sqrt_r, evaluatedConvo, evaluatedPowU0);
    fptype Ipy   = a0 * Ipy_0 + a1 * Ipy_1 + a2 * Ipy_2 + a3 * Ipy_3;
    fptype Imy   = a0 * Imy_0 + a1 * Imy_1 + a2 * Imy_2 + a3 * Imy_3;
    // fudge factor 10 by comparison with old values. Need to really understand what is missing
    _P1 += 1. * preFactor * (preFactor_plus * Ipy + preFactor_minus * Imy);
    _P3 += 1. * preFactor * (preFactor_plus * Ipy - preFactor_minus * Imy);
}

__device__ fptype EvaluateAcceptanceGn(int n, int k, fptype _sqrt_r, fptype *_evaluatedConvo, fptype *_evaluatedPowU) {
    // n : I_gn
    // k: order of polynmoial term
    fptype res = 0.;
    for(int i = 0; i <= n + k; i++) {
        res += BinomCoeff(n + k, i) * _evaluatedConvo[n + k - i] * _evaluatedPowU[i];
    }
    // is this term double-counted? ALready considered in EvaluateKnotSinCos? -> no, it's correct
    res = res * pow(1. / _sqrt_r, k);
    return res;
}

__device__ void EvaluateKnotSinCos(fptype &_P2,
                                   fptype &_P4,
                                   fptype Gamma,
                                   fptype knot_low,
                                   fptype knot_high,
                                   fptype a0,
                                   fptype a1,
                                   fptype a2,
                                   fptype a3,
                                   fptype tprime,
                                   fptype sigma,
                                   fptype x) {
    // sqrt_r = 1/sigma
    fptype sqrt_r = 1. / sigma;
    // r = 1/sigma^2
    fptype r  = sqrt_r * sqrt_r;
    fptype p  = 2 * Gamma - 2 * (tprime)*r;
    fptype q  = (tprime) * (tprime)*r;
    fptype u0 = p / (2. * sqrt_r);
    // factor 1/sqrt(2PI) absorbed in small phi?
    fptype preFactor    = 1. / sqrt_r / sigma * exp(-0.5 * (q - u0 * u0));
    fptype commonFactor = (x * Gamma / sqrt_r);

    fptype evaluatedConvo[9];
    fptype evaluatedPowU0[9];
    fptype lowLim  = sqrt_r * knot_low + u0;
    fptype highLim = sqrt_r * knot_high + u0;
    for(int i = 0; i < 9; i++) {
        evaluatedConvo[i] = EvaluateConvo(i, lowLim, highLim);
    }
    evaluatedPowU0[0] = 1.;
    for(int i = 1; i < 9; i++) {
        evaluatedPowU0[i] = evaluatedPowU0[i - 1] * -u0;
    }
    fptype Ig0_0 = EvaluateAcceptanceGn(0, 0, sqrt_r, evaluatedConvo, evaluatedPowU0);
    fptype Ig0_1 = EvaluateAcceptanceGn(0, 1, sqrt_r, evaluatedConvo, evaluatedPowU0);
    fptype Ig0_2 = EvaluateAcceptanceGn(0, 2, sqrt_r, evaluatedConvo, evaluatedPowU0);
    fptype Ig0_3 = EvaluateAcceptanceGn(0, 3, sqrt_r, evaluatedConvo, evaluatedPowU0);

    fptype Ig1_0 = EvaluateAcceptanceGn(1, 0, sqrt_r, evaluatedConvo, evaluatedPowU0);
    fptype Ig1_1 = EvaluateAcceptanceGn(1, 1, sqrt_r, evaluatedConvo, evaluatedPowU0);
    fptype Ig1_2 = EvaluateAcceptanceGn(1, 2, sqrt_r, evaluatedConvo, evaluatedPowU0);
    fptype Ig1_3 = EvaluateAcceptanceGn(1, 3, sqrt_r, evaluatedConvo, evaluatedPowU0);

    fptype Ig2_0 = EvaluateAcceptanceGn(2, 0, sqrt_r, evaluatedConvo, evaluatedPowU0);
    fptype Ig2_1 = EvaluateAcceptanceGn(2, 1, sqrt_r, evaluatedConvo, evaluatedPowU0);
    fptype Ig2_2 = EvaluateAcceptanceGn(2, 2, sqrt_r, evaluatedConvo, evaluatedPowU0);
    fptype Ig2_3 = EvaluateAcceptanceGn(2, 3, sqrt_r, evaluatedConvo, evaluatedPowU0);

    fptype Ig3_0 = EvaluateAcceptanceGn(3, 0, sqrt_r, evaluatedConvo, evaluatedPowU0);
    fptype Ig3_1 = EvaluateAcceptanceGn(3, 1, sqrt_r, evaluatedConvo, evaluatedPowU0);
    fptype Ig3_2 = EvaluateAcceptanceGn(3, 2, sqrt_r, evaluatedConvo, evaluatedPowU0);
    fptype Ig3_3 = EvaluateAcceptanceGn(3, 3, sqrt_r, evaluatedConvo, evaluatedPowU0);

    fptype Ig4_0 = EvaluateAcceptanceGn(4, 0, sqrt_r, evaluatedConvo, evaluatedPowU0);
    fptype Ig4_1 = EvaluateAcceptanceGn(4, 1, sqrt_r, evaluatedConvo, evaluatedPowU0);
    fptype Ig4_2 = EvaluateAcceptanceGn(4, 2, sqrt_r, evaluatedConvo, evaluatedPowU0);
    fptype Ig4_3 = EvaluateAcceptanceGn(4, 3, sqrt_r, evaluatedConvo, evaluatedPowU0);

    fptype Ig5_0 = EvaluateAcceptanceGn(5, 0, sqrt_r, evaluatedConvo, evaluatedPowU0);
    fptype Ig5_1 = EvaluateAcceptanceGn(5, 1, sqrt_r, evaluatedConvo, evaluatedPowU0);
    fptype Ig5_2 = EvaluateAcceptanceGn(5, 2, sqrt_r, evaluatedConvo, evaluatedPowU0);
    fptype Ig5_3 = EvaluateAcceptanceGn(5, 3, sqrt_r, evaluatedConvo, evaluatedPowU0);

    fptype Ig0 = a0 * Ig0_0 + a1 * Ig0_1 + a2 * Ig0_2 + a3 * Ig0_3;
    // Ig0 *= pow(commonFactor, 0);
    fptype Ig1 = a0 * Ig1_0 + a1 * Ig1_1 + a2 * Ig1_2 + a3 * Ig1_3;
    Ig1 *= commonFactor;
    fptype Ig2 = a0 * Ig2_0 + a1 * Ig2_1 + a2 * Ig2_2 + a3 * Ig2_3;
    commonFactor *= commonFactor;
    Ig2 *= commonFactor;
    fptype Ig3 = a0 * Ig3_0 + a1 * Ig3_1 + a2 * Ig3_2 + a3 * Ig3_3;
    commonFactor *= commonFactor;
    Ig3 *= commonFactor;
    fptype Ig4 = a0 * Ig4_0 + a1 * Ig4_1 + a2 * Ig4_2 + a3 * Ig4_3;
    commonFactor *= commonFactor;
    Ig4 *= commonFactor;
    fptype Ig5 = a0 * Ig5_0 + a1 * Ig5_1 + a2 * Ig5_2 + a3 * Ig5_3;
    commonFactor *= commonFactor;
    Ig5 *= commonFactor;

    // fudge factor 20 by comparison with old values. Need to really understand what is missing -> had a factor 0.5 in
    // preFactor which should not have been there
    _P2 += preFactor * (Ig0 - 0.5 * Ig2 + 1. / 24. * Ig4);
    _P4 += preFactor * (Ig1 - 1. / 6. * Ig3 + 1. / 120. * Ig5);
}

__device__ void gaussian_splice(fptype &_P1,
                                fptype &_P2,
                                fptype &_P3,
                                fptype &_P4,
                                fptype _tau,
                                fptype adjTime,
                                fptype xmixing,
                                fptype ymixing,
                                fptype adjSigma,
                                int nKnots,
                                fptype *knots,
                                fptype *spline_0,
                                fptype *spline_1,
                                fptype *spline_2,
                                fptype *spline_3) {
    fptype Gamma = 1. / _tau;
    _P1          = 0.;
    _P2          = 0.;
    _P3          = 0.;
    _P4          = 0.;
    for(int i = 0; i < nKnots - 1; i++) {
        EvaluateKnot(_P1,
                     _P3,
                     Gamma,
                     knots[i],
                     knots[i + 1],
                     spline_0[i],
                     spline_1[i],
                     spline_2[i],
                     spline_3[i],
                     adjTime,
                     adjSigma,
                     ymixing);
        EvaluateKnotSinCos(_P2,
                           _P4,
                           Gamma,
                           knots[i],
                           knots[i + 1],
                           spline_0[i],
                           spline_1[i],
                           spline_2[i],
                           spline_3[i],
                           adjTime,
                           adjSigma,
                           xmixing);
    }
}

__device__ void mytmpgaussian(fptype &_P1,
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

__device__ void mygaussian_high(fptype &_P1,
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
    fptype _1oSigma  = 1. / adjSigma;
    fptype _1o2SqrtA = 0.5 * _1oSqrtA; // 1/(2 sqrt(r))
    fptype _1oSigma2 = _1oSigma * _1oSigma;
    fptype _NormG    = SQRT1o2PI * _1oSigma;

    fptype _q   = 0.5 * adjTime * adjTime * _1oSigma2;
    fptype _Bgn = -adjTime * _1oSigma2;

    fptype _Gamma = 1 / _tau;
    fptype _B     = _Gamma + _Bgn;

    fptype _u0  = _1o2SqrtA * (_B + selbias);
    fptype u1   = 1. / _1oSqrtA * Tthres + _u0;
    fptype _u02 = _u0 * _u0;
    fptype u12  = u1 * u1;
    fptype _F   = _1oSqrtA * exp(-_q + _u02);

    fptype _Ig0 = SQRTPIo2 * erfc(u1);
    fptype _Ig1 = 0.5 * exp(-1 * u1 * u1);
    fptype _Ig2 = _Ig1 * u1 + 0.5 * _Ig0;
    fptype _Ig3 = _Ig1 * (u12 + 1);

    fptype _R   = xmixing * _Gamma * _1oSqrtA;
    fptype _R2  = _R * _R;
    fptype _It0 = _F * _Ig0;
    fptype _It1 = _F * _R * (_Ig1 - _u0 * _Ig0);
    fptype _It2 = _F * _R2 * (_Ig2 - _u0 * _Ig1 * 2 + _u02 * _Ig0);
    fptype _It3 = _F * _R2 * _R * (_Ig3 - _u0 * _Ig2 * 3 + _u02 * _Ig1 * 3 - _u0 * _u02 * _Ig0);

    // fptype _P0   = _NormG *  _It0;
    _P2 += preConst * _NormG * (_It0 - 0.5 * _It2);
    _P4 += preConst * _NormG * (_It1 - R1o6 * _It3);

    fptype _u0py = _1o2SqrtA * (_B - ymixing * _Gamma + selbias);
    fptype _u0my = _1o2SqrtA * (_B + ymixing * _Gamma + selbias);
    fptype u1py  = _u0py + 1. / _1oSqrtA * Tthres;
    fptype u1my  = _u0my + 1. / _1oSqrtA * Tthres;
    fptype _Fpy  = _1oSqrtA * exp(-_q + _u0py * _u0py);
    fptype _Fmy  = _1oSqrtA * exp(-_q + _u0my * _u0my);
    fptype _Ipy  = _Fpy * SQRTPIo2 * erfc(u1py);
    fptype _Imy  = _Fmy * SQRTPIo2 * erfc(u1my);
    _P1 += preConst * _NormG * 0.5 * (_Ipy + _Imy);
    _P3 += preConst * _NormG * 0.5 * (_Ipy - _Imy);
}

__device__ fptype device_threegauss_resolutionSplice(fptype coshterm,
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

    int nKnots   = pc.getConstant(0);
    int nSplines = nKnots - 1;
    // TODO: upper splice end should probably not be to large, as function might become negative
    // TODO: at the moment support up to 8 knots/ 7 splines!
    fptype knots[8];
    fptype spline_0[8];
    fptype spline_1[7];
    fptype spline_2[7];
    fptype spline_3[7];
    for(int i = 0; i < nKnots; i++) {
        knots[i] = pc.getParameter(8 + i);
    }
    knots[0] = 0.;
    for(int i = 0; i < nSplines; i++) {
        int offset  = 8 + nKnots + i;
        spline_0[i] = pc.getParameter(offset);
        spline_1[i] = pc.getParameter(offset + nSplines);
        spline_2[i] = pc.getParameter(offset + 2 * nSplines);
        spline_3[i] = pc.getParameter(offset + 3 * nSplines);
    }
    // find exp with same slope
    fptype m        = spline_1[nSplines - 1];
    fptype b        = spline_0[nSplines - 1];
    fptype lastKnot = knots[nKnots - 1];
    fptype fVal     = m * lastKnot + b;
    fptype fDeriv   = m;
    fptype expSlope = -fDeriv / fVal;
    fptype expConst = fVal;

    fptype cp1;
    fptype cp2;
    fptype cp3;
    fptype cp4;

    fptype tp1;
    fptype tp2;
    fptype tp3;
    fptype tp4;

    fptype op1;
    fptype op2;
    fptype op3;
    fptype op4;

    gaussian_splice(cp1,
                    cp2,
                    cp3,
                    cp4,
                    tau,
                    dtime - coreBias * sigma,
                    xmixing,
                    ymixing,
                    coreScaleFactor * sigma,
                    nKnots,
                    knots,
                    spline_0,
                    spline_1,
                    spline_2,
                    spline_3);
    gaussian_splice(tp1,
                    tp2,
                    tp3,
                    tp4,
                    tau,
                    dtime - tailBias * sigma,
                    xmixing,
                    ymixing,
                    tailScaleFactor * sigma,
                    nKnots,
                    knots,
                    spline_0,
                    spline_1,
                    spline_2,
                    spline_3);
    gaussian_splice(op1,
                    op2,
                    op3,
                    op4,
                    tau,
                    dtime - outlBias * sigma,
                    xmixing,
                    ymixing,
                    outlScaleFactor * sigma,
                    nKnots,
                    knots,
                    spline_0,
                    spline_1,
                    spline_2,
                    spline_3);

    mygaussian_high(cp1,
                    cp2,
                    cp3,
                    cp4,
                    tau,
                    dtime - coreBias * sigma,
                    xmixing,
                    ymixing,
                    coreScaleFactor * sigma,
                    expSlope,
                    lastKnot,
                    expConst);
    mygaussian_high(tp1,
                    tp2,
                    tp3,
                    tp4,
                    tau,
                    dtime - tailBias * sigma,
                    xmixing,
                    ymixing,
                    tailScaleFactor * sigma,
                    expSlope,
                    lastKnot,
                    expConst);
    mygaussian_high(op1,
                    op2,
                    op3,
                    op4,
                    tau,
                    dtime - outlBias * sigma,
                    xmixing,
                    ymixing,
                    outlScaleFactor * sigma,
                    expSlope,
                    lastKnot,
                    expConst);

    /*
    printf("------------ \n");
    printf("CP1: %f %f %f \n", cp1_tmp, cp1, cp1/cp1_tmp);
    printf("TP1: %f %f %f \n", tp1_tmp, op1, tp1/tp1_tmp);
    printf("OP1: %f %f %f \n", op1_tmp, tp1, op1/op1_tmp);

    printf("CP2: %f %f %f \n", cp2_tmp, cp2, cp2/cp2_tmp);
    printf("TP2: %f %f %f \n", tp2_tmp, op2, tp2/tp2_tmp);
    printf("OP2: %f %f %f \n", op2_tmp, tp2, op2/op2_tmp);

    printf("CP3: %f %f %f \n", cp3_tmp, cp3, cp3/cp3_tmp);
    printf("TP3: %f %f %f \n", tp3_tmp, op3, tp3/tp3_tmp);
    printf("OP3: %f %f %f \n", op3_tmp, tp3, op3/op3_tmp);

    printf("CP4: %f %f %f \n", cp4_tmp, cp4, cp4/cp4_tmp);
    printf("TP4: %f %f %f \n", tp4_tmp, op4, tp4/tp4_tmp);
    printf("OP4: %f %f %f \n", op4_tmp, tp4, op4/op4_tmp);
    printf("------------ \n");
    */

    fptype _P1 = coreFraction * (cp1) + tailFraction * (tp1) + outlFraction * (op1);
    fptype _P2 = coreFraction * (cp2) + tailFraction * (tp2) + outlFraction * (op2);
    fptype _P3 = coreFraction * (cp3) + tailFraction * (tp3) + outlFraction * (op3);
    fptype _P4 = coreFraction * (cp4) + tailFraction * (tp4) + outlFraction * (op4);

    fptype ret = 0;
    ret += coshterm * _P1;
    ret += costerm * _P2;
    ret -= 2 * sinhterm * _P3;
    ret -= 2 * sinterm * _P4; // Notice sign difference wrt to Mikhail's code, because I have AB* and he has A*B.

    // pc.incrementIndex (1, 8, 0, 0, 1);

    // cuPrintf("device_threegauss_resolution %f %f %f %f %f\n", coshterm, costerm, sinhterm, sinterm, dtime);
    // printf("Splcie ret: %f %f %f %f %f \n", ret, _P1, _P2, _P3, _P4);
    return ret;
}

__device__ device_resfunction_ptr ptr_to_threegaussSplice = device_threegauss_resolutionSplice;

ThreeGaussResolutionSplice::ThreeGaussResolutionSplice(Variable cf,
                                                       Variable tf,
                                                       Variable cb,
                                                       Variable cs,
                                                       Variable tb,
                                                       Variable ts,
                                                       Variable ob,
                                                       Variable os,
                                                       std::vector<Variable> knots,
                                                       std::vector<Variable> a0,
                                                       std::vector<Variable> a1,
                                                       std::vector<Variable> a2,
                                                       std::vector<Variable> a3)
    : MixingTimeResolution("ThreeGaussResolutionSplice", cf, tf, cb, cs, tb, ts, ob, os)
    , m_knots(knots)
    , m_a0(a0)
    , m_a1(a1)
    , m_a2(a2)
    , m_a3(a3) {
    initIndex();
    registerConstant(knots.size());
    for(auto knot : knots)
        registerParameter(knot);
    for(auto i : a0)
        registerParameter(i);
    for(auto i : a1)
        registerParameter(i);
    for(auto i : a2)
        registerParameter(i);
    for(auto i : a3)
        registerParameter(i);

    registerFunction("ptr_to_threegaussSplice", ptr_to_threegaussSplice);
}
ThreeGaussResolutionSplice::~ThreeGaussResolutionSplice() = default;

__host__ __device__ fptype
NormTermY(fptype tlow, fptype thigh, fptype yterm, fptype Gamma, int k, fptype *evaluatedPowGammaYTerm) {
    // k: order of polynomial term
    fptype res = 0.;
    for(int i = 0; i <= k; i++) {
        fptype preFactor = -1. / evaluatedPowGammaYTerm[i + 1] * Factorial(k) / Factorial(k - i);
        fptype curIntegral
            = (pow(thigh, k - i) * exp(-thigh * Gamma * yterm)) - (pow(tlow, k - i) * exp(-tlow * Gamma * yterm));
        res += preFactor * curIntegral;
    }
    return res;
}

/*
fptype NormTermX(fptype tlow, fptype thigh, fpcomplex xterm, fptype Gamma, int k) {
    //k: order of polynomial term
    fpcomplex res{0.,0.};
    for (int i = 0; i <= k; i++) {
        fptype preFactor = -1./pow(Gamma*xterm, i+1) * Factorial(k)/Factorial(k-i);
        fptype curIntegral = (pow(thigh, k-i) * exp(-thigh * Gamma * xterm)) - (pow(tlow, k-i) * exp(-tlow * Gamma *
xterm)); res += preFactor * curIntegral;
    }

    return res.real();
}
*/

fptype ThreeGaussResolutionSplice::normalization(
    fptype di1, fptype di2, fptype di3, fptype di4, fptype tau, fptype xmixing, fptype ymixing) const {
    fptype Gamma = 1. / tau;
    // NB! In thesis notation, A_1 = (A + B), A_2 = (A - B).
    // Here di1 = |A^2|, di2 = |B^2|, di3,4 = Re,Im(AB^*).
    // Distinction between numerical subscribts and A,B is crucial
    // for comparing thesis math to this math!

    // fptype timeIntegralOne = tau / (1 - ymixing * ymixing);
    // fptype timeIntegralTwo = tau / (1 + xmixing * xmixing);
    // fptype timeIntegralThr = ymixing * timeIntegralOne;
    // fptype timeIntegralFou = xmixing * timeIntegralTwo;

    auto nKnots = m_knots.size();
    // upper splice end should probably not be to large, as function might become negative
    std::vector<fptype> knots;
    std::vector<fptype> spline_0;
    std::vector<fptype> spline_1;
    std::vector<fptype> spline_2;
    std::vector<fptype> spline_3;
    for(auto i : m_knots) {
        knots.push_back(i.getValue());
    }
    for(auto i : m_a0)
        spline_0.push_back(i.getValue());
    for(auto i : m_a1)
        spline_1.push_back(i.getValue());
    for(auto i : m_a2)
        spline_2.push_back(i.getValue());
    for(auto i : m_a3)
        spline_3.push_back(i.getValue());

    /*copying the code from ThreeGaussResolutionExt
    the exponential there is of the form expCosnt * exp(-expSlope(t - Tthres))
    */
    auto nSplines   = nKnots - 1;
    fptype m        = spline_1[nSplines - 1];
    fptype b        = spline_0[nSplines - 1];
    fptype lastKnot = knots[nKnots - 1];
    fptype fVal     = m * lastKnot + b;
    fptype fDeriv   = m;
    fptype expSlope = -fDeriv / fVal;
    // can simplify
    // fptype expConst = fVal/exp(-expSlope * lastKnot)/exp(expSlope * lastKnot);
    fptype expConst = fVal;

    /*
    0., 0.2,  0.35, 0.5,  0.65, 1.05, 2.1,  4.5
    fptype  spline_0[] = {0.03765277000655794, 0.03765277000655794, 0.033653622225172486, 0.05183326055919843,
    0.037589800853886655, 0.03943035933470114, 0.03461351985201154}; fptype  spline_1[] = {0.008125428187878276,
    0.00812542818787828, 0.04240383774261078, -0.06667399226154494, -0.0009349474677982864, -0.006193685984411092,
    -0.005360727973783324}; fptype  spline_2[] = {-0.026572145610579714, -0.02657214561057973, -0.12451045862410115,
    0.09364520138421031, -0.007491790606169145, -0.002483468209395046, 0.0}; fptype  spline_3[] =
    {0.0027387047655195807, 0.0027387047655195998, 0.09601328858792094, -0.049423818084286704, 0.0024413060133437765,
    0.000851362395320254, 0.0};
    */

    fptype timeIntegralOne   = 0.;
    fptype timeIntegralTwo   = 0.;
    fptype timeIntegralThree = 0.;
    fptype timeIntegralFour  = 0.;
    fptype evaluatedPowGammaYPlus[10];
    fptype evaluatedPowGammaYMinus[10];
    fptype evaluatedPowGamma[10];
    evaluatedPowGammaYPlus[0]  = 1.;
    evaluatedPowGammaYMinus[0] = 1.;
    evaluatedPowGamma[0]       = 1.;
    for(int i = 1; i < 10; i++) {
        evaluatedPowGammaYPlus[i]  = evaluatedPowGammaYPlus[i - 1] * Gamma * (ymixing + 1.);
        evaluatedPowGammaYMinus[i] = evaluatedPowGammaYMinus[i - 1] * Gamma * (1. - ymixing);
        evaluatedPowGamma[i]       = evaluatedPowGamma[i - 1] * Gamma;
    }

    // 1. - ymixing or ymixing - 1.
    for(int i = 0; i < nKnots - 1; i++) {
        timeIntegralOne
            += spline_0.at(i)
               * (0.5 * NormTermY(knots.at(i), knots.at(i + 1), 1. - ymixing, Gamma, 0, evaluatedPowGammaYMinus)
                  + 0.5 * NormTermY(knots.at(i), knots.at(i + 1), ymixing + 1., Gamma, 0, evaluatedPowGammaYPlus));
        timeIntegralOne
            += spline_1.at(i)
               * (0.5 * NormTermY(knots.at(i), knots.at(i + 1), 1. - ymixing, Gamma, 1, evaluatedPowGammaYMinus)
                  + 0.5 * NormTermY(knots.at(i), knots.at(i + 1), ymixing + 1., Gamma, 1, evaluatedPowGammaYPlus));
        timeIntegralOne
            += spline_2.at(i)
               * (0.5 * NormTermY(knots.at(i), knots.at(i + 1), 1. - ymixing, Gamma, 2, evaluatedPowGammaYMinus)
                  + 0.5 * NormTermY(knots.at(i), knots.at(i + 1), ymixing + 1., Gamma, 2, evaluatedPowGammaYPlus));
        timeIntegralOne
            += spline_3.at(i)
               * (0.5 * NormTermY(knots.at(i), knots.at(i + 1), 1. - ymixing, Gamma, 3, evaluatedPowGammaYMinus)
                  + 0.5 * NormTermY(knots.at(i), knots.at(i + 1), ymixing + 1., Gamma, 3, evaluatedPowGammaYPlus));

        timeIntegralThree
            += spline_0.at(i)
               * (0.5 * NormTermY(knots.at(i), knots.at(i + 1), 1. - ymixing, Gamma, 0, evaluatedPowGammaYMinus)
                  - 0.5 * NormTermY(knots.at(i), knots.at(i + 1), ymixing + 1., Gamma, 0, evaluatedPowGammaYPlus));
        timeIntegralThree
            += spline_1.at(i)
               * (0.5 * NormTermY(knots.at(i), knots.at(i + 1), 1. - ymixing, Gamma, 1, evaluatedPowGammaYMinus)
                  - 0.5 * NormTermY(knots.at(i), knots.at(i + 1), ymixing + 1., Gamma, 1, evaluatedPowGammaYPlus));
        timeIntegralThree
            += spline_2.at(i)
               * (0.5 * NormTermY(knots.at(i), knots.at(i + 1), 1. - ymixing, Gamma, 2, evaluatedPowGammaYMinus)
                  - 0.5 * NormTermY(knots.at(i), knots.at(i + 1), ymixing + 1., Gamma, 2, evaluatedPowGammaYPlus));
        timeIntegralThree
            += spline_3.at(i)
               * (0.5 * NormTermY(knots.at(i), knots.at(i + 1), 1. - ymixing, Gamma, 3, evaluatedPowGammaYMinus)
                  - 0.5 * NormTermY(knots.at(i), knots.at(i + 1), ymixing + 1., Gamma, 3, evaluatedPowGammaYPlus));

        timeIntegralTwo += spline_0.at(i)
                           * (NormTermY(knots.at(i), knots.at(i + 1), 1., Gamma, 0, evaluatedPowGamma)
                              - 0.5 * pow(xmixing * Gamma, 2)
                                    * NormTermY(knots.at(i), knots.at(i + 1), 1., Gamma, 2, evaluatedPowGamma)
                              + 1. / 24. * pow(xmixing * Gamma, 4)
                                    * NormTermY(knots.at(i), knots.at(i + 1), 1., Gamma, 4, evaluatedPowGamma));
        timeIntegralTwo += spline_1.at(i)
                           * (NormTermY(knots.at(i), knots.at(i + 1), 1., Gamma, 1, evaluatedPowGamma)
                              - 0.5 * pow(xmixing * Gamma, 2)
                                    * NormTermY(knots.at(i), knots.at(i + 1), 1., Gamma, 3, evaluatedPowGamma)
                              + 1. / 24. * pow(xmixing * Gamma, 4)
                                    * NormTermY(knots.at(i), knots.at(i + 1), 1., Gamma, 5, evaluatedPowGamma));
        timeIntegralTwo += spline_2.at(i)
                           * (NormTermY(knots.at(i), knots.at(i + 1), 1., Gamma, 2, evaluatedPowGamma)
                              - 0.5 * pow(xmixing * Gamma, 2)
                                    * NormTermY(knots.at(i), knots.at(i + 1), 1., Gamma, 4, evaluatedPowGamma)
                              + 1. / 24. * pow(xmixing * Gamma, 4)
                                    * NormTermY(knots.at(i), knots.at(i + 1), 1., Gamma, 6, evaluatedPowGamma));
        timeIntegralTwo += spline_3.at(i)
                           * (NormTermY(knots.at(i), knots.at(i + 1), 1., Gamma, 3, evaluatedPowGamma)
                              - 0.5 * pow(xmixing * Gamma, 2)
                                    * NormTermY(knots.at(i), knots.at(i + 1), 1., Gamma, 5, evaluatedPowGamma)
                              + 1. / 24. * pow(xmixing * Gamma, 4)
                                    * NormTermY(knots.at(i), knots.at(i + 1), 1., Gamma, 7, evaluatedPowGamma));

        timeIntegralFour
            += spline_0.at(i)
               * (xmixing * Gamma * NormTermY(knots.at(i), knots.at(i + 1), 1., Gamma, 1, evaluatedPowGamma)
                  - 1. / 6. * pow(xmixing * Gamma, 3)
                        * NormTermY(knots.at(i), knots.at(i + 1), 1., Gamma, 3, evaluatedPowGamma)
                  + 1. / 120. * pow(xmixing * Gamma, 5)
                        * NormTermY(knots.at(i), knots.at(i + 1), 1., Gamma, 5, evaluatedPowGamma));
        timeIntegralFour
            += spline_1.at(i)
               * (xmixing * Gamma * NormTermY(knots.at(i), knots.at(i + 1), 1., Gamma, 2, evaluatedPowGamma)
                  - 1. / 6. * pow(xmixing * Gamma, 3)
                        * NormTermY(knots.at(i), knots.at(i + 1), 1., Gamma, 4, evaluatedPowGamma)
                  + 1. / 120. * pow(xmixing * Gamma, 5)
                        * NormTermY(knots.at(i), knots.at(i + 1), 1., Gamma, 6, evaluatedPowGamma));
        timeIntegralFour
            += spline_2.at(i)
               * (xmixing * Gamma * NormTermY(knots.at(i), knots.at(i + 1), 1., Gamma, 3, evaluatedPowGamma)
                  - 1. / 6. * pow(xmixing * Gamma, 3)
                        * NormTermY(knots.at(i), knots.at(i + 1), 1., Gamma, 5, evaluatedPowGamma)
                  + 1. / 120. * pow(xmixing * Gamma, 5)
                        * NormTermY(knots.at(i), knots.at(i + 1), 1., Gamma, 7, evaluatedPowGamma));
        timeIntegralFour
            += spline_3.at(i)
               * (xmixing * Gamma * NormTermY(knots.at(i), knots.at(i + 1), 1., Gamma, 4, evaluatedPowGamma)
                  - 1. / 6. * pow(xmixing * Gamma, 3)
                        * NormTermY(knots.at(i), knots.at(i + 1), 1., Gamma, 6, evaluatedPowGamma)
                  + 1. / 120. * pow(xmixing * Gamma, 5)
                        * NormTermY(knots.at(i), knots.at(i + 1), 1., Gamma, 8, evaluatedPowGamma));

        /*
        fpcomplex imag_x = {0., xmixing}
         for (int i = 0; i < nKnots-1; i++ ) {
        }
        timeIntegralTwo += spline_0.at(i) * (0.5 * NormTermX(knots.at(i), knots.at(i+1), 1. - imag_x, Gamma, 0) + 0.5 *
        NormTermX(knots.at(i), knots.at(i+1), imag_x + 1., Gamma, 0) ); timeIntegralTwo += spline_1.at(i) * (0.5 *
        NormTermX(knots.at(i), knots.at(i+1), 1. - imag_x, Gamma, 1) + 0.5 * NormTermX(knots.at(i), knots.at(i+1),
        imag_x + 1., Gamma, 1) ); timeIntegralTwo += spline_2.at(i) * (0.5 * NormTermX(knots.at(i), knots.at(i+1), 1. -
        imag_x, Gamma, 2) + 0.5 * NormTermX(knots.at(i), knots.at(i+1), imag_x+ 1., Gamma, 2) ); timeIntegralTwo +=
        spline_3.at(i) * (0.5 * NormTermX(knots.at(i), knots.at(i+1), 1. - imag_x, Gamma, 3) + 0.5 *
        NormTermX(knots.at(i), knots.at(i+1), imag_x+ 1., Gamma, 3) );

        timeIntegralFour += spline_0.at(i) * (0.5 * NormTermX(knots.at(i), knots.at(i+1), 1. - imag_x, Gamma, 0) - 0.5 *
        NormTermX(knots.at(i), knots.at(i+1), imag_x + 1., Gamma, 0) ); timeIntegralFour += spline_1.at(i) * (0.5 *
        NormTermX(knots.at(i), knots.at(i+1), 1. - imag_x, Gamma, 1) - 0.5 * NormTermX(knots.at(i), knots.at(i+1),
        imag_x + 1., Gamma, 1) ); timeIntegralFour += spline_2.at(i) * (0.5 * NormTermX(knots.at(i), knots.at(i+1), 1. -
        imag_x, Gamma, 2) - 0.5 * NormTermX(knots.at(i), knots.at(i+1), imag_x + 1., Gamma, 2) ); timeIntegralFour +=
        spline_3.at(i) * (0.5 * NormTermX(knots.at(i), knots.at(i+1), 1. - imag_x, Gamma, 3) - 0.5 *
        NormTermX(knots.at(i), knots.at(i+1), imag_x + 1., Gamma, 3) );
        */
    }

    fptype selBias_high  = expSlope;
    fptype Tthres        = knots[nKnots - 1];
    fptype preConst_high = expConst * exp(selBias_high * Tthres);

    fptype gammaPlusBias = Gamma + selBias_high;

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

    timeIntegralOne += timeIntegralOne_high;
    timeIntegralTwo += timeIntegralTwo_high;
    timeIntegralThree += timeIntegralThr_high;
    timeIntegralFour += timeIntegralFour_high;

    fptype ret = timeIntegralOne * (di1 + di2); // ~ |A|^2 + |B|^2
    ret += timeIntegralTwo * (di1 - di2);       // ~ |A|^2 - |B|^2
    ret -= 2 * timeIntegralThree * di3;         // ~ Re(A_1 A_2^*)
    ret -= 2 * timeIntegralFour * di4;          // ~ Im(A_1 A_2^*)
    // printf("Splcie norm: %f \n", ret);

    /*
    fptype selBias = 0.21868332187637923;
    fptype timeIntegralOne_old
        = (selBias + 1 / tau) / (selBias * selBias + 2 * selBias / tau + (1 - ymixing * ymixing) / (tau * tau));
    fptype timeIntegralTwo_old
        = (selBias + 1 / tau) / (selBias * selBias + 2 * selBias / tau + (1 + xmixing * xmixing) / (tau * tau));
    fptype timeIntegralThr_old
        = (ymixing / tau) / (selBias * selBias + 2 * selBias / tau + (1 - ymixing * ymixing) / (tau * tau));
    fptype timeIntegralFou_old
        = (xmixing / tau) / (selBias * selBias + 2 * selBias / tau + (1 + xmixing * xmixing) / (tau * tau));



    printf("------------ \n");
    printf("timeIntegralOne: %f %f %f \n", timeIntegralOne_old, timeIntegralOne, timeIntegralOne/timeIntegralOne_old);
    printf("timeIntegralTwo: %f %f %f \n", timeIntegralTwo_old, timeIntegralTwo, timeIntegralTwo/timeIntegralTwo_old);
    printf("timeIntegralThree: %f %f %f \n", timeIntegralThr_old, timeIntegralThree,
    timeIntegralThree/timeIntegralThr_old); printf("timeIntegralFour: %f %f %f \n", timeIntegralFou_old,
    timeIntegralFour, timeIntegralFour/timeIntegralFou_old); printf("------------ \n");
     */
    return ret;
}

} // namespace GooFit
