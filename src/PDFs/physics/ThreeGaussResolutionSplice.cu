#include <cmath>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/ThreeGaussResolutionSplice.h>

namespace GooFit {

const fptype R1o6 = 1.0 / 6.0;
#define SQRTPIo2 (1.0 / M_2_SQRTPI)
#define SQRT1o2PI (sqrt(0.5 * M_1_PI))
#define SQRT1oPI (sqrt(M_1_PI))



__device__ fptype BigPhi (fptype u0, fptype u1) {
  fptype res = 0.5*(erfc(u0 * SQRT1oPI) - erfc(u1 * SQRT1oPI));
  return res;
}

__device__ fptype SmallPhi(fptype u0, fptype u1) {
    fptype res = SQRT1o2PI * (exp(-0.5*u1*u1) - exp(-0.5*u0*u0));
    return res;
}

__device__ int DoubleFactorial(int j) {
  int res = 1.;
  while (j > 0) {
    res = res * j;
    j = j -2;
  }
  return res;
}

__device__ fptype FactorialOdd(int n, fptype u0, fptype u1) {
    int k = (n - 1)/2;
    int constFactorial = DoubleFactorial(2*k);
    fptype res = 0.;
    for (int j = 0; j <= k; j++) {
        res = res + 1./DoubleFactorial(j) * (pow(u1, 2*j) - pow(u0, 2*j));
    }
    res = res * -SmallPhi(u0, u1) * constFactorial;
    return res;
}

__device__ fptype FactorialEven(int n, fptype u0, fptype u1) {
    int k = (n - 2)/2;
    int constFactorial = DoubleFactorial(2*k + 1);
    fptype res = 0.;
    for (int j = 0; j <= k; j++) {
        res = res + 1./DoubleFactorial(2*j + 1) * (pow(u1, 2*j + 1) - pow(u0, 2*j + 1));
    }
   res = res * -SmallPhi(u0, u1) * constFactorial;
   res = res + DoubleFactorial(k*2 + 1)*BigPhi(u0, u1);
   return res;
}

__device__ fptype EvaluateConvo(int n, fptype u0, fptype u1) {
  if(n==0) return BigPhi(u0, u1);
  if(n==1) return -SmallPhi(u0, u1);
  if(n%2 == 0) return FactorialEven(n, u0, u1);
  else return FactorialOdd(n, u0, u1);
}

__device__ int Factorial(int j) {
   if (j == 0) return 1;
   int res = 1.;
   while(j > 0) {
    res = res * j;
    j--;
   }
   return j;
}

__device__ int BinomCoeff(int n, int k) {
    int res = Factorial(n);
    res = res /(Factorial(k) * Factorial(n-k));
    return res;
}

__device__ fptype EvaluateAcceptance(int n, fptype r, fptype u0, fptype tlow, fptype thigh) {
    fptype res = 0.;
    for(int k = 0; k <=n; k++) {
      res += BinomCoeff(n, k) * EvaluateConvo(n-k, sqrt(r)*tlow + u0, sqrt(r)*thigh + u0) * pow(-u0,k);
    }
    res = res * pow(1./sqrt(r), n);
    return res;
}

__device__ void EvaluateKnot(fptype &_P1, fptype &_P3, fptype Gamma, fptype knot_low, fptype knot_high, fptype a0, fptype a1, fptype a2, fptype a3, fptype tprime, fptype sigma, fptype y) {
    fptype mu = 0.;
    fptype r = 1./sigma/sigma;
    fptype p_plus  = 2*Gamma*(1.-y) - (tprime - mu)/sigma/sigma;
    fptype p_minus = 2*Gamma*(1.+y) - (tprime - mu)/sigma/sigma;
    fptype q = (tprime - mu)/sigma/sigma;
    fptype u0_plus  = p_plus/(2.*sqrt(r));
    fptype u0_minus = p_minus/(2.*sqrt(r));

    fptype preFactor = 0.5/sqrt(2.*M_PI*r)/sigma;
    fptype preFactor_plus = exp(-0.5*(q - u0_plus*u0_plus));
    fptype preFactor_minus = exp(-0.5*(q - u0_minus*u0_minus));
    fptype Ipy_0  = EvaluateAcceptance(0, r, u0_plus, knot_low, knot_high);
    fptype Ipy_1  = EvaluateAcceptance(1, r, u0_plus, knot_low, knot_high);
    fptype Ipy_2  = EvaluateAcceptance(2, r, u0_plus, knot_low, knot_high);
    fptype Ipy_3  = EvaluateAcceptance(3, r, u0_plus, knot_low, knot_high);

    fptype Imy_0  = EvaluateAcceptance(0, r, u0_minus, knot_low, knot_high);
    fptype Imy_1  = EvaluateAcceptance(1, r, u0_minus, knot_low, knot_high);
    fptype Imy_2  = EvaluateAcceptance(2, r, u0_minus, knot_low, knot_high);
    fptype Imy_3  = EvaluateAcceptance(3, r, u0_minus, knot_low, knot_high);
    fptype Ipy = a0 * Ipy_0 + a1 * Ipy_1  + a2 * Ipy_2 + a3 * Ipy_3 ;
    fptype Imy = a0 * Imy_0 + a1 * Imy_1  + a2 * Imy_2 + a3 * Imy_3 ;
    _P1 += preFactor * (preFactor_plus * Ipy + preFactor_minus * Imy);
    _P3 += preFactor * (preFactor_plus * Ipy - preFactor_minus * Imy);

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
                         fptype * knots,
                         fptype * spline_0,
                         fptype * spline_1,
                         fptype * spline_2,
                         fptype * spline_3) {
    fptype Gamma = 1./_tau;
    _P1 = 0.;
    _P2 = 0.;
    _P3 = 0.;
    _P4 = 0.;
    for (int i = 0; i < nKnots-1; i++) {
       EvaluateKnot(_P1, _P3, Gamma, knots[i], knots[i+1], spline_0[i], spline_1[i], spline_2[i], spline_3[i], adjTime, adjSigma, ymixing); 
    }
    
    
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
    fptype tailFraction    = (1 - coreFraction) * pc.getParameter(1);
    fptype outlFraction    = 1 - coreFraction - tailFraction;
    fptype coreBias        = pc.getParameter(2);
    fptype coreScaleFactor = pc.getParameter(3);
    fptype tailBias        = pc.getParameter(4);
    fptype tailScaleFactor = pc.getParameter(5);
    fptype outlBias        = pc.getParameter(6);
    fptype outlScaleFactor = pc.getParameter(7);
    fptype selbias_low         = pc.getParameter(8);
    fptype selbias_high         = pc.getParameter(9);
    fptype Tthreshold  = pc.getParameter(10);
    fptype constantC = pc.getParameter(11);

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
    int nKnots = 8;
    fptype knots[] = {0.05, 0.2,  0.35, 0.5,  0.65, 1.05, 2.1,  2.7};
    fptype  spline_0[] = {3.86541588E-2,  3.89678917E-2,  3.81065346E-2,  3.64586010E-2,
   3.51188071E-2,  3.18690217E-2,  2.41612462E-2};
    fptype  spline_1[] = {6.24580670E-3, -1.94405803E-3, -9.42205861E-3, -1.05535214E-2,
  -8.11840539E-3, -8.15136576E-3, -5.59415052E-3};
    fptype  spline_2[] = {-2.84859894E-2, -2.61131088E-2, -2.37402283E-2,  1.61971427E-2,
   3.69643757E-5, -1.19365309E-4,  0.0};
    fptype  spline_3[] = {5.27306786E-3,  5.27306786E-3,  8.87497133E-2, -3.59115074E-2,
  -1.30274737E-4,  8.48944035E-4,  0.0};
      

    gaussian_splice(cp1_low, cp2_low, cp3_low, cp4_low, tau, dtime - coreBias * sigma, xmixing, ymixing, coreScaleFactor * sigma, nKnots, knots, spline_0, spline_1, spline_2, spline_3);
    gaussian_splice(tp1_low, tp2_low, tp3_low, tp4_low, tau, dtime - tailBias * sigma, xmixing, ymixing, tailScaleFactor * sigma, nKnots, knots, spline_0, spline_1, spline_2, spline_3);
    gaussian_splice(op1_low, op2_low, op3_low, op4_low, tau, dtime - outlBias * sigma, xmixing, ymixing, outlScaleFactor * sigma, nKnots, knots, spline_0, spline_1, spline_2, spline_3);




    fptype _P1 = coreFraction * (cp1_low) + tailFraction * (tp1_low) + outlFraction * (op1_low);
    fptype _P2 = coreFraction * (cp2_low) + tailFraction * (tp2_low) + outlFraction * (op2_low);
    fptype _P3 = coreFraction * (cp3_low) + tailFraction * (tp3_low) + outlFraction * (op3_low);
    fptype _P4 = coreFraction * (cp4_low) + tailFraction * (tp4_low) + outlFraction * (op4_low);

    fptype ret = 0;
    ret += coshterm * _P1;
    ret += costerm * _P2;
    ret -= 2 * sinhterm * _P3;
    ret -= 2 * sinterm * _P4; // Notice sign difference wrt to Mikhail's code, because I have AB* and he has A*B.

    // pc.incrementIndex (1, 8, 0, 0, 1);

    // cuPrintf("device_threegauss_resolution %f %f %f %f %f\n", coshterm, costerm, sinhterm, sinterm, dtime);
    return ret;
}

__device__ device_resfunction_ptr ptr_to_threegaussSplice = device_threegauss_resolutionSplice;

ThreeGaussResolutionSplice::ThreeGaussResolutionSplice(
    Variable cf, Variable tf, Variable cb, Variable cs, Variable tb, Variable ts, Variable ob, Variable os, Variable sb_low, Variable sb_high, Variable Tthres, Variable constantC)
    : MixingTimeResolution("ThreeGaussResolutionSplice", cf, tf, cb, cs, tb, ts, ob, os, sb_low, sb_high, Tthres, constantC)
    , selectionBias_low(sb_low), selectionBias_high(sb_high), mTthreshold(Tthres), mConstantC(constantC) {
    initIndex();

    registerFunction("ptr_to_threegaussSplice", ptr_to_threegaussSplice);
}
ThreeGaussResolutionSplice::~ThreeGaussResolutionSplice() = default;

fptype ThreeGaussResolutionSplice::normalization(
    fptype di1, fptype di2, fptype di3, fptype di4, fptype tau, fptype xmixing, fptype ymixing) const {
    
    // NB! In thesis notation, A_1 = (A + B), A_2 = (A - B).
    // Here di1 = |A^2|, di2 = |B^2|, di3,4 = Re,Im(AB^*).
    // Distinction between numerical subscribts and A,B is crucial
    // for comparing thesis math to this math!

    // fptype timeIntegralOne = tau / (1 - ymixing * ymixing);
    // fptype timeIntegralTwo = tau / (1 + xmixing * xmixing);
    // fptype timeIntegralThr = ymixing * timeIntegralOne;
    // fptype timeIntegralFou = xmixing * timeIntegralTwo;

    fptype selBias_low = selectionBias_low.getValue();
    fptype selBias_high = selectionBias_high.getValue();
    fptype Tthres = mTthreshold.getValue();
    fptype C = mConstantC.getValue();
    fptype preConst_low = C * exp(selBias_low * Tthres);
    fptype preConst_high = C * exp(selBias_high * Tthres);

    fptype Gamma = 1./tau; 
    fptype gammaPlusBias = Gamma + selBias_low;

    fptype timeIntegralOne_low = 0.5*preConst_low * (  1./(ymixing*Gamma - Gamma - selBias_low)  * (   exp( Tthres * (ymixing*Gamma - Gamma - selBias_low) )  - 1  )  +
        1./(-ymixing*Gamma - Gamma -selBias_low)  * (   exp( Tthres * (-ymixing*Gamma - Gamma - selBias_low) )  - 1  )  );

    fptype timeIntegralThr_low = 0.5*preConst_low * (  1./(ymixing*Gamma - Gamma - selBias_low)  * (   exp( Tthres * (ymixing*Gamma - Gamma - selBias_low) )  - 1  )  -
        1./(-ymixing*Gamma - Gamma -selBias_low)  * (   exp( Tthres * (-ymixing*Gamma - Gamma - selBias_low) )  - 1  )  );

    fptype timeIntegralTwo_low = preConst_low * (exp(-gammaPlusBias * Tthres) / (gammaPlusBias*gammaPlusBias + xmixing*xmixing *Gamma*Gamma) * ( -gammaPlusBias * cos(xmixing *Gamma * Tthres)  + xmixing*Gamma*sin(xmixing*Gamma*Tthres)) - (-gammaPlusBias)/(gammaPlusBias*gammaPlusBias + xmixing*xmixing *Gamma*Gamma) );

    fptype timeIntegralFour_low = preConst_low * (exp(-gammaPlusBias * Tthres) / (gammaPlusBias*gammaPlusBias + xmixing*xmixing *Gamma*Gamma) * ( -gammaPlusBias * sin(xmixing *Gamma * Tthres)  - xmixing*Gamma*cos(xmixing*Gamma*Tthres)) - (-xmixing*Gamma)/(gammaPlusBias*gammaPlusBias + xmixing*xmixing *Gamma*Gamma) );


    gammaPlusBias = Gamma + selBias_high;

    fptype timeIntegralOne_high = 0.5*preConst_high * (  -1./(ymixing*Gamma - Gamma - selBias_high)  * (   exp( Tthres * (ymixing*Gamma - Gamma - selBias_high) )   )  -
        1./(-ymixing*Gamma - Gamma -selBias_high)  * (   exp( Tthres * (-ymixing*Gamma - Gamma - selBias_high) )   )  );

    fptype timeIntegralThr_high = 0.5*preConst_high * (  -1./(ymixing*Gamma - Gamma - selBias_high)  * (   exp( Tthres * (ymixing*Gamma - Gamma - selBias_high) )    )  +
        1./(-ymixing*Gamma - Gamma -selBias_high)  * (   exp( Tthres * (-ymixing*Gamma - Gamma - selBias_high) )   )  );


    fptype timeIntegralTwo_high = preConst_high * (-exp(-gammaPlusBias * Tthres) / (gammaPlusBias*gammaPlusBias + xmixing*xmixing *Gamma*Gamma) * ( -gammaPlusBias * cos(xmixing *Gamma * Tthres)  + xmixing*Gamma*sin(xmixing*Gamma*Tthres))  );

    fptype timeIntegralFour_high = preConst_high * (-exp(-gammaPlusBias * Tthres) / (gammaPlusBias*gammaPlusBias + xmixing*xmixing *Gamma*Gamma) * ( -gammaPlusBias * sin(xmixing *Gamma * Tthres)  - xmixing*Gamma*cos(xmixing*Gamma*Tthres))  );




    fptype timeIntegralOne =  timeIntegralOne_low + timeIntegralOne_high;
    fptype timeIntegralTwo = timeIntegralTwo_low + timeIntegralTwo_high;
    fptype timeIntegralThr = timeIntegralThr_low + timeIntegralThr_high;
    fptype timeIntegralFou = timeIntegralFour_low + timeIntegralFour_high;

    fptype ret = timeIntegralOne * (di1 + di2); // ~ |A|^2 + |B|^2
    ret += timeIntegralTwo * (di1 - di2);        // ~ |A|^2 - |B|^2
    ret -= 2 * timeIntegralThr * di3;          // ~ Re(A_1 A_2^*)
    ret -= 2 * timeIntegralFou * di4;           // ~ Im(A_1 A_2^*)

    return ret;
}
} // namespace GooFit
