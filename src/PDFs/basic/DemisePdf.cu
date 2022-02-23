#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/basic/DemisePdf.h>
#include <goofit/Variable.h>
#include <stdio.h>

namespace GooFit {

__device__ auto device_Demise(fptype *evt, ParameterContainer &pc) -> fptype {
    int id = pc.getObservable(0);

    fptype x     = RO_CACHE(evt[id]);
    fptype m0    = pc.getParameter(0);
    fptype sigma = pc.getParameter(1);

    pc.incrementIndex(1, 2, 0, 1, 1);

    double t = x - m0;
    fptype arg = t/sigma;
    
    fptype ret = (t >= 0)?0:1.-exp(arg);

    return ret;
}


__device__ device_function_ptr ptr_to_Demise = device_Demise;


__host__ DemisePdf::DemisePdf(std::string n, Observable _x, Variable m0, Variable sigma)
    : GooPdf("DemisePdf", n, _x, m0, sigma) {

//  This is a PDF that falls exponentially to zero as the argument approaches its
//  maximum value (m0) but is almost  one at low values of the argument.
//
//  Below the maximum value it has the form [1 - exp((x-m0)/sigma)]

    registerFunction("ptr_to_Demise", ptr_to_Demise);

    initialize();
}

__host__ auto DemisePdf::integrate(fptype lo, fptype hi) const -> fptype {
    fptype m0    = host_parameters[parametersIdx + 1];
    fptype sigma = host_parameters[parametersIdx + 2];

    fptype norm = 0;

    if (hi > m0) {
     norm = (hi-lo) - sigma;}
    else {
     norm = (hi-lo) - sigma*exp((hi-m0)/(lo-m0));
    }

    return norm;
}

} // namespace GooFit
