#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/basic/CrystalBallPdf.h>
#include <goofit/Variable.h>

namespace GooFit {

__device__ fptype device_CrystalBall(fptype *evt, ParameterContainer &pc) {
    // Left-hand tail if alpha is less than 0,
    // right-hand tail if greater, pure Gaussian if 0.
    // return 1;
    int id = pc.getObservable(0);

    fptype x     = RO_CACHE(evt[id]);
    fptype mean  = pc.getParameter(0);
    fptype sigma = pc.getParameter(1);
    fptype alpha = pc.getParameter(2);
    fptype power = pc.getParameter(3);
    fptype rx    = (sigma != 0) ? (x - mean) / sigma : 0;
    fptype ret   = 0;

    if((alpha > 0 && rx <= alpha) || // Right-hand tail, in Gaussian region
       (alpha < 0 && rx >= alpha) || // Left-hand tail, in Gaussian region
       (alpha == 0)) {               // Pure Gaussian
        ret = exp(-0.5 * rx * rx);
    } else { // Tail part
        fptype n_over_alpha = power / alpha;
        fptype a            = exp(-0.5 * alpha * alpha);
        fptype b            = n_over_alpha - alpha;
        fptype d            = b + rx;
        d                   = (d != 0) ? n_over_alpha / d : 0;
        ret                 = a * pow(d, power);
    }

    pc.incrementIndex(1, 4, 0, 1, 1);

    // if ((0 == THREADIDX) && (0 == BLOCKIDX)) printf("device_CB: %f %f %f %f %f %f\n", x, mean, sigma, alpha, power,
    // ret);
    return ret;
}

__device__ device_function_ptr ptr_to_CrystalBall = device_CrystalBall;

__host__ CrystalBallPdf::CrystalBallPdf(std::string n, Observable _x, Variable mean, Variable sigma, Variable alpha)
    : CrystalBallPdf(n, _x, mean, sigma, alpha, Variable(n + "_n", 2)) {}

__host__ CrystalBallPdf::CrystalBallPdf(
    std::string n, Observable _x, Variable mean, Variable sigma, Variable alpha, Variable power)
    : GooPdf("CrystalBallPdf", n, _x, mean, sigma, alpha, power) {
    registerFunction("ptr_to_CrystalBall", ptr_to_CrystalBall);

    initialize();
}

__host__ fptype CrystalBallPdf::integrate(fptype lo, fptype hi) const {
    static const fptype sqrtPiOver2 = 1.2533141373;
    static const fptype sqrt2       = 1.4142135624;

    fptype result = 0.0;
    bool useLog   = false;

    fptype mean  = parametersList[0].getValue();
    fptype sigma = parametersList[1].getValue();
    fptype alpha = parametersList[2].getValue();
    fptype power = parametersList[3].getValue();

    if(fabs(power - 1.0) < 1.0e-05)
        useLog = true;

    fptype tmin = (lo - mean) / sigma;
    fptype tmax = (hi - mean) / sigma;

    if(alpha < 0) {
        fptype tmp = tmin;
        tmin       = -tmax;
        tmax       = -tmp;
    }

    fptype absAlpha = fabs(alpha);

    if(tmin >= -absAlpha) {
        result += sigma * sqrtPiOver2 * (erf(tmax / sqrt2) - erf(tmin / sqrt2));
    } else if(tmax <= -absAlpha) {
        fptype a = pow(power / absAlpha, power) * exp(-0.5 * absAlpha * absAlpha);
        fptype b = power / absAlpha - absAlpha;

        if(useLog) {
            result += a * sigma * (log(b - tmin) - log(b - tmax));
        } else {
            result += a * sigma / (1.0 - power)
                      * (1.0 / (pow(b - tmin, power - 1.0)) - 1.0 / (pow(b - tmax, power - 1.0)));
        }
    } else {
        fptype a = pow(power / absAlpha, power) * exp(-0.5 * absAlpha * absAlpha);
        fptype b = power / absAlpha - absAlpha;

        fptype term1 = 0.0;

        if(useLog) {
            term1 = a * sigma * (log(b - tmin) - log(power / absAlpha));
        } else {
            term1 = a * sigma / (1.0 - power)
                    * (1.0 / (pow(b - tmin, power - 1.0)) - 1.0 / (pow(power / absAlpha, power - 1.0)));
        }

        fptype term2 = sigma * sqrtPiOver2 * (erf(tmax / sqrt2) - erf(-absAlpha / sqrt2));
        result += term1 + term2;
    }

    return result;
}

} // namespace GooFit
