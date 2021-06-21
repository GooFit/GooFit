#include <goofit/Log.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/basic/BernsteinPdf.h>
#include <goofit/Variable.h>

namespace GooFit {

__device__ auto factorial(int n) -> int {
    if(n == 0)
        return 1;
    else
        return n * factorial(n - 1);
}

__device__ auto binomial(int k, int n) -> fptype {
    if(k > n)
        return 0;
    else {
        fptype a = factorial(n);
        fptype b = factorial(k) * factorial(n - k);
        return a / b;
    }
}

//__device__ fptype device_BernsteinPdf (fptype* evt, fptype* p, unsigned int* indices) {
__device__ auto device_BernsteinPdf(fptype *evt, ParameterContainer &pc) -> fptype {
    // nP (lowX uppX ord) c0 .. cn nO (o1 o2 o3 o4)

    int id = pc.getObservable(0);

    int numParams     = pc.getNumParameters();
    int num_constants = pc.getNumConstants();

    fptype x = RO_CACHE(evt[id]);

    fptype lowX         = pc.getConstant(1);
    fptype uppX         = pc.getConstant(2);
    unsigned int mindeg = pc.getConstant(3);
    int maxDegree       = pc.getConstant(0) + mindeg;
    // unsigned int maxDegreePlus = maxDegree;//;+1;

    fptype xNorm = (x - lowX) / (uppX - lowX);

    fptype ret = 0;

    for(int i = 0; i < maxDegree - mindeg; i++) {
        // printf("Eval: %d - %d \n",i,maxDegree);
        fptype val = pc.getParameter(i);
        // fptype par = val;
        fptype bin = binomial(i, maxDegree - 1);

        val *= bin * pow(xNorm, i) * pow(1 - xNorm, maxDegree - 1 - i);
        val /= (uppX - lowX);
        ret += val;

        // printf("x = %f part = %f val = %f xNorm = %f i = %d bin = %f ret = %f \n",x,par,val,xNorm,i,bin,ret);
    }

    // printf("ret = %f \n",ret);
    pc.incrementIndex(1, numParams, num_constants, 1, 1);
    return ret;
}

__device__ device_function_ptr ptr_to_BernsteinPdf = device_BernsteinPdf;

__host__ BernsteinPdf::BernsteinPdf(std::string n, Observable obs, std::vector<Variable> coeffs, unsigned int mindeg)
    : GooPdf("BernsteinPdf", n, obs) {
    // Only 1-dimensional!

    // unsigned int numParameters = (maxDegree + 1);

    registerConstant(int(coeffs.size()));
    registerConstant(obs.getLowerLimit());
    registerConstant(obs.getUpperLimit());
    registerConstant(mindeg);

    for(auto &coeff : coeffs) {
        registerParameter(coeff);
        // printf("%s %f %f %f \n",coeff.getName(),coeff.getValue(),coeff.getLowerLimit(),coeff.getUpperLimit());
        // printf("%d\n",i);
        // i++;
    }

    registerFunction("ptr_to_BernsteinPdf", ptr_to_BernsteinPdf);

    initialize();
}

} // namespace GooFit
