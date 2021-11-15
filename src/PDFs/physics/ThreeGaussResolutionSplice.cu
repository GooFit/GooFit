#include <cmath>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/ThreeGaussResolutionSplice.h>

namespace GooFit {

const fptype R1o6 = 1.0 / 6.0;
#define SQRTPIo2 (1.0 / M_2_SQRTPI)
#define SQRT1o2PI (sqrt(0.5 * M_1_PI))
#define SQRT1oPI (sqrt(M_1_PI))



__device__ fptype BigPhi (fptype u0, fptype u1) {
  fptype res = 0.5*(erfc(u0 /sqrt(2.)) - erfc(u1 /sqrt(2.)));
  return res;
}

__device__ fptype SmallPhi(fptype u0, fptype u1) {
    fptype res = SQRT1o2PI * (exp(-0.5*u1*u1) - exp(-0.5*u0*u0));
    return res;
}


__device__ fptype SmallPhi(fptype x) {
    fptype res = SQRT1o2PI * (exp(-0.5*x*x));
    return res;
}


__device__ int DoubleFactorial(int j) {
  //printf("DoubleFactorial start \n");
  int res = 1.;
  while (j > 0) {
    res = res * j;
    j = j -2;
  }
  //printf("DoubleFactorial finish \n");
  return res;
}

__device__ fptype FactorialOdd(int n, fptype u0, fptype u1) {
    //evaluate odd moments of Gaussian
   // printf("FactorialOdd start \n");
    int k = (n - 1)/2;
    int constFactorial = DoubleFactorial(2*k);
    fptype res = 0.;
    for (int j = 0; j <= k; j++) {
        res = res + 1./DoubleFactorial(j) * (-SmallPhi(u1)*pow(u1, 2*j) + SmallPhi(u0)*pow(u0, 2*j));
    }
    res = res * constFactorial;
    //printf("FactorialOdd finish \n");
    return res;
}

__device__ fptype FactorialEven(int n, fptype u0, fptype u1) {
    //evaluate even moments of Gaussian
   // printf("FactorialOdd start \n");
    int k = (n - 2)/2;
    int constFactorial = DoubleFactorial(2*k + 1);
    fptype res = 0.;
    for (int j = 0; j <= k; j++) {
        res = res + 1./DoubleFactorial(2*j + 1) * (-SmallPhi(u1)*pow(u1, 2*j + 1) + SmallPhi(u0)*pow(u0, 2*j + 1));
    }
   res = res *  constFactorial;
   res = res + constFactorial*BigPhi(u0, u1);
   //printf("FactorialOdd finish \n");
   return res;
}

__device__ fptype EvaluateConvo(int n, fptype u0, fptype u1) {
  //evaluate moments of Gaussian function
  if(n==0) return BigPhi(u0, u1);
  if(n==1) return -SmallPhi(u0, u1);
  if(n%2 == 0) return FactorialEven(n, u0, u1);
  else return FactorialOdd(n, u0, u1);
}

__device__ int Factorial(int j) {
   if (j == 0) return 1;
  // printf("evaluating factorial of %d \n", j);
   int res = 1.;
   while(j > 0) {
    res = res * j;
   // printf("factorial res is %d \n", res);
    j--;
   }
   return res;
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
    fptype p_plus  = 2*Gamma*(1.-y) - 2*(tprime - mu)/sigma/sigma;
    fptype p_minus = 2*Gamma*(1.+y) - 2*(tprime - mu)/sigma/sigma;
    fptype q = (tprime - mu)*(tprime - mu)/sigma/sigma;
    fptype u0_plus  = p_plus/(2.*sqrt(r));
    fptype u0_minus = p_minus/(2.*sqrt(r));
    //factor 1/sqrt(2PI) absorbed in small phi?
    fptype preFactor = 0.5/sqrt(r)/sigma;
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


__device__ fptype EvaluateAcceptanceGn(int n, int k, fptype r, fptype u0, fptype tlow, fptype thigh) {
    //n : I_gn
    //k: order of polynmoial term
    fptype res = 0.;
    for(int i = 0; i <=n+k; i++) {
      res += BinomCoeff(n+k, i) * EvaluateConvo(n+k-i, sqrt(r)*tlow + u0, sqrt(r)*thigh + u0) * pow(-u0,i);
    }
    res = res * pow(1./sqrt(r), k);
    return res;
}


__device__ void EvaluateKnotSinCos(fptype &_P2, fptype &_P4, fptype Gamma, fptype knot_low, fptype knot_high, fptype a0, fptype a1, fptype a2, fptype a3, fptype tprime, fptype sigma, fptype x) {
    fptype mu = 0.;
    fptype r = 1./sigma/sigma;
    fptype p  = 2*Gamma - 2*(tprime - mu)/sigma/sigma;
    fptype q = (tprime - mu)*(tprime - mu)/sigma/sigma;
    fptype u0  = p/(2.*sqrt(r));
    //factor 1/sqrt(2PI) absorbed in small phi?
    fptype preFactor = 0.5/sqrt(r)/sigma * exp(-0.5*(q - u0*u0));
    fptype commonFactor = (x * Gamma /sqrt(r));
    fptype Ig0_0  = EvaluateAcceptanceGn(0, 0, r, u0, knot_low, knot_high);
    fptype Ig0_1  = EvaluateAcceptanceGn(0, 1, r, u0, knot_low, knot_high);
    fptype Ig0_2  = EvaluateAcceptanceGn(0, 1, r, u0, knot_low, knot_high);
    fptype Ig0_3  = EvaluateAcceptanceGn(0, 3, r, u0, knot_low, knot_high);

    fptype Ig1_0  = EvaluateAcceptanceGn(1, 0, r, u0, knot_low, knot_high);
    fptype Ig1_1  = EvaluateAcceptanceGn(1, 1, r, u0, knot_low, knot_high);
    fptype Ig1_2  = EvaluateAcceptanceGn(1, 1, r, u0, knot_low, knot_high);
    fptype Ig1_3  = EvaluateAcceptanceGn(1, 3, r, u0, knot_low, knot_high);

    fptype Ig2_0  = EvaluateAcceptanceGn(2, 0, r, u0, knot_low, knot_high);
    fptype Ig2_1  = EvaluateAcceptanceGn(2, 1, r, u0, knot_low, knot_high);
    fptype Ig2_2  = EvaluateAcceptanceGn(2, 1, r, u0, knot_low, knot_high);
    fptype Ig2_3  = EvaluateAcceptanceGn(2, 3, r, u0, knot_low, knot_high);

    fptype Ig3_0  = EvaluateAcceptanceGn(3, 0, r, u0, knot_low, knot_high);
    fptype Ig3_1  = EvaluateAcceptanceGn(3, 1, r, u0, knot_low, knot_high);
    fptype Ig3_2  = EvaluateAcceptanceGn(3, 1, r, u0, knot_low, knot_high);
    fptype Ig3_3  = EvaluateAcceptanceGn(3, 3, r, u0, knot_low, knot_high);


    
    fptype Ig0 = a0 * Ig0_0 + a1 * Ig0_1  + a2 * Ig0_2 + a3 * Ig0_3;
    //Ig0 *= pow(commonFactor, 0);
    fptype Ig1 = a0 * Ig1_0 + a1 * Ig1_1  + a2 * Ig1_2 + a3 * Ig1_3;
    Ig1 *= commonFactor;
    fptype Ig2 = a0 * Ig2_0 + a1 * Ig2_1  + a2 * Ig2_2 + a3 * Ig2_3;
    Ig2 *= pow(commonFactor, 2);
    fptype Ig3 = a0 * Ig3_0 + a1 * Ig3_1  + a2 * Ig3_2 + a3 * Ig3_3;
    Ig3 *= pow(commonFactor, 3);

    _P2 += preFactor * (Ig0 - 0.5 * Ig2);
    _P4 += preFactor * (Ig1 - 1./6. * Ig3);

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
       EvaluateKnotSinCos(_P2, _P4, Gamma, knots[i], knots[i+1], spline_0[i], spline_1[i], spline_2[i], spline_3[i], adjTime, adjSigma, xmixing);
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
 



    int nKnots = pc.getConstant(0);
    int nSplines = nKnots - 1;
    // TODO: upper splice end should probably not be to large, as fucntion might become negative
    //TODO: at the moment support up to 8 knots/ 7 splines!
    fptype knots[8];
    fptype spline_0[8];
    fptype spline_1[7];
    fptype spline_2[7];
    fptype spline_3[7];
    for(int i = 0; i < nKnots; i++) {
        knots[i] = pc.getParameter(8+i);

    }
    for(int i = 0; i < nSplines; i++) {
        int offset = 8 + nKnots + i;
        spline_0[i] = pc.getParameter(offset);
        spline_1[i] = pc.getParameter(offset + nSplines);
        spline_2[i] = pc.getParameter(offset + 2*nSplines);
        spline_3[i] = pc.getParameter(offset + 3*nSplines);

    }


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

    gaussian_splice(cp1, cp2, cp3, cp4, tau, dtime - coreBias * sigma, xmixing, ymixing, coreScaleFactor * sigma, nKnots, knots, spline_0, spline_1, spline_2, spline_3);
    gaussian_splice(tp1, tp2, tp3, tp4, tau, dtime - tailBias * sigma, xmixing, ymixing, tailScaleFactor * sigma, nKnots, knots, spline_0, spline_1, spline_2, spline_3);
    gaussian_splice(op1, op2, op3, op4, tau, dtime - outlBias * sigma, xmixing, ymixing, outlScaleFactor * sigma, nKnots, knots, spline_0, spline_1, spline_2, spline_3);




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
    //printf("Splcie ret: %f %f %f %f %f \n", ret, _P1, _P2, _P3, _P4);
    return ret;
}

__device__ device_resfunction_ptr ptr_to_threegaussSplice = device_threegauss_resolutionSplice;

ThreeGaussResolutionSplice::ThreeGaussResolutionSplice(
    Variable cf, Variable tf, Variable cb, Variable cs, Variable tb, Variable ts, Variable ob, Variable os, std::vector<Variable> knots, std::vector<Variable> a0, std::vector<Variable> a1, std::vector<Variable> a2, std::vector<Variable> a3)
    : MixingTimeResolution("ThreeGaussResolutionSplice", cf, tf, cb, cs, tb, ts, ob, os)
    , m_knots(knots), m_a0(a0), m_a1(a1), m_a2(a2), m_a3(a3) {
    initIndex();
    registerConstant(knots.size());
    for (auto knot : knots) registerParameter(knot);
    for (auto i : a0) registerParameter(i);
    for (auto i : a1) registerParameter(i);
    for (auto i : a2) registerParameter(i);
    for (auto i : a3) registerParameter(i);

    registerFunction("ptr_to_threegaussSplice", ptr_to_threegaussSplice);
}
ThreeGaussResolutionSplice::~ThreeGaussResolutionSplice() = default;






fptype NormTermY(fptype tlow, fptype thigh, fptype yterm, fptype Gamma, int k) {
    //k: order of polynomial term
    fptype res = 0.;
    for (int i = 0; i <= k; i++) {
        fptype preFactor = -1./pow(Gamma*yterm, i+1) * Factorial(k)/Factorial(k-i);
        fptype curIntegral = (pow(thigh, k-i) * exp(-thigh * Gamma * yterm)) - (pow(tlow, k-i) * exp(-tlow * Gamma * yterm));
        res += preFactor * curIntegral;
    }
    return res;
}


fptype ThreeGaussResolutionSplice::normalization(
    fptype di1, fptype di2, fptype di3, fptype di4, fptype tau, fptype xmixing, fptype ymixing) const {

    fptype Gamma = 1./tau;
    // NB! In thesis notation, A_1 = (A + B), A_2 = (A - B).
    // Here di1 = |A^2|, di2 = |B^2|, di3,4 = Re,Im(AB^*).
    // Distinction between numerical subscribts and A,B is crucial
    // for comparing thesis math to this math!

    // fptype timeIntegralOne = tau / (1 - ymixing * ymixing);
    // fptype timeIntegralTwo = tau / (1 + xmixing * xmixing);
    // fptype timeIntegralThr = ymixing * timeIntegralOne;
    // fptype timeIntegralFou = xmixing * timeIntegralTwo;

    auto nKnots = m_knots.size();
    // upper splice end should probably not be to large, as fucntion might become negative
    std::vector<fptype> knots;
    std::vector<fptype>  spline_0;
    std::vector<fptype>  spline_1;
    std::vector<fptype>  spline_2;
    std::vector<fptype>  spline_3;
    for(auto i : m_knots) {
        knots.push_back(i.getValue());
    }
    for(auto i : m_a0) spline_0.push_back(i.getValue());
    for(auto i : m_a1) spline_1.push_back(i.getValue());
    for(auto i : m_a2) spline_2.push_back(i.getValue());
    for(auto i : m_a3) spline_3.push_back(i.getValue());

    /*
    0., 0.2,  0.35, 0.5,  0.65, 1.05, 2.1,  4.5
    fptype  spline_0[] = {0.03765277000655794, 0.03765277000655794, 0.033653622225172486, 0.05183326055919843, 0.037589800853886655, 0.03943035933470114, 0.03461351985201154};
    fptype  spline_1[] = {0.008125428187878276, 0.00812542818787828, 0.04240383774261078, -0.06667399226154494, -0.0009349474677982864, -0.006193685984411092, -0.005360727973783324};
    fptype  spline_2[] = {-0.026572145610579714, -0.02657214561057973, -0.12451045862410115, 0.09364520138421031, -0.007491790606169145, -0.002483468209395046, 0.0};
    fptype  spline_3[] = {0.0027387047655195807, 0.0027387047655195998, 0.09601328858792094, -0.049423818084286704, 0.0024413060133437765, 0.000851362395320254, 0.0};
    */
    





    fptype timeIntegralOne =  0.;
    fptype timeIntegralTwo =  0.;
    fptype timeIntegralThree =  0.;
    fptype timeIntegralFour =  0.;
    // 1. - ymixing or ymixing - 1.
    for (int i = 0; i < nKnots-1; i++ ) {
        timeIntegralOne += spline_0.at(i) * (0.5 * NormTermY(knots.at(i), knots.at(i+1), 1. - ymixing, Gamma, 0) + 0.5 * NormTermY(knots.at(i), knots.at(i+1), ymixing + 1., Gamma, 0) );
        timeIntegralOne += spline_1.at(i) * (0.5 * NormTermY(knots.at(i), knots.at(i+1), 1. - ymixing, Gamma, 1) + 0.5 * NormTermY(knots.at(i), knots.at(i+1), ymixing + 1., Gamma, 1) );
        timeIntegralOne += spline_2.at(i) * (0.5 * NormTermY(knots.at(i), knots.at(i+1), 1. - ymixing, Gamma, 2) + 0.5 * NormTermY(knots.at(i), knots.at(i+1), ymixing + 1., Gamma, 2) );
        timeIntegralOne += spline_3.at(i) * (0.5 * NormTermY(knots.at(i), knots.at(i+1), 1. - ymixing, Gamma, 3) + 0.5 * NormTermY(knots.at(i), knots.at(i+1), ymixing + 1., Gamma, 3) );

        timeIntegralThree += spline_0.at(i) * (0.5 * NormTermY(knots.at(i), knots.at(i+1), 1. - ymixing, Gamma, 0) - 0.5 * NormTermY(knots.at(i), knots.at(i+1), ymixing + 1., Gamma, 0) );
        timeIntegralThree += spline_1.at(i) * (0.5 * NormTermY(knots.at(i), knots.at(i+1), 1. - ymixing, Gamma, 1) - 0.5 * NormTermY(knots.at(i), knots.at(i+1), ymixing + 1., Gamma, 1) );
        timeIntegralThree += spline_2.at(i) * (0.5 * NormTermY(knots.at(i), knots.at(i+1), 1. - ymixing, Gamma, 2) - 0.5 * NormTermY(knots.at(i), knots.at(i+1), ymixing + 1., Gamma, 2) );
        timeIntegralThree += spline_3.at(i) * (0.5 * NormTermY(knots.at(i), knots.at(i+1), 1. - ymixing, Gamma, 3) - 0.5 * NormTermY(knots.at(i), knots.at(i+1), ymixing + 1., Gamma, 3) );

        timeIntegralTwo += spline_0.at(i) * (NormTermY(knots.at(i), knots.at(i+1), 1., Gamma, 0) - 0.5 * pow(xmixing*Gamma,2) * NormTermY(knots.at(i), knots.at(i+1), 1., Gamma, 2) );
        timeIntegralTwo += spline_1.at(i) * (NormTermY(knots.at(i), knots.at(i+1), 1., Gamma, 1) - 0.5 * pow(xmixing*Gamma,2) * NormTermY(knots.at(i), knots.at(i+1), 1., Gamma, 3) );
        timeIntegralTwo += spline_2.at(i) * (NormTermY(knots.at(i), knots.at(i+1), 1., Gamma, 2) - 0.5 * pow(xmixing*Gamma,2) * NormTermY(knots.at(i), knots.at(i+1), 1., Gamma, 4) );
        timeIntegralTwo += spline_3.at(i) * (NormTermY(knots.at(i), knots.at(i+1), 1., Gamma, 3) - 0.5 * pow(xmixing*Gamma,2) * NormTermY(knots.at(i), knots.at(i+1), 1., Gamma, 5) );

        timeIntegralFour += spline_0.at(i) * (xmixing*Gamma*NormTermY(knots.at(i), knots.at(i+1), 1., Gamma, 1) - 1./6. * pow(xmixing*Gamma,3) * NormTermY(knots.at(i), knots.at(i+1), 1., Gamma, 3) );
        timeIntegralFour += spline_1.at(i) * (xmixing*Gamma*NormTermY(knots.at(i), knots.at(i+1), 1., Gamma, 2) - 1./6. * pow(xmixing*Gamma,3) * NormTermY(knots.at(i), knots.at(i+1), 1., Gamma, 4) );
        timeIntegralFour += spline_2.at(i) * (xmixing*Gamma*NormTermY(knots.at(i), knots.at(i+1), 1., Gamma, 3) - 1./6. * pow(xmixing*Gamma,3) * NormTermY(knots.at(i), knots.at(i+1), 1., Gamma, 5) );
        timeIntegralFour += spline_3.at(i) * (xmixing*Gamma*NormTermY(knots.at(i), knots.at(i+1), 1., Gamma, 4) - 1./6. * pow(xmixing*Gamma,3) * NormTermY(knots.at(i), knots.at(i+1), 1., Gamma, 6) );

    }

    fptype ret = timeIntegralOne * (di1 + di2); // ~ |A|^2 + |B|^2
    ret += timeIntegralTwo * (di1 - di2);        // ~ |A|^2 - |B|^2
    ret -= 2 * timeIntegralThree * di3;          // ~ Re(A_1 A_2^*)
    ret -= 2 * timeIntegralFour * di4;           // ~ Im(A_1 A_2^*)
    //printf("Splcie norm: %f \n", ret);
    return ret;


}


} // namespace GooFit
