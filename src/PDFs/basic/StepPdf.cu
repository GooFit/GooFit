#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/basic/StepPdf.h>
#include <iostream>

//  modified 210923 Mike Sokoloff
//  added argument Up to signature; if its value  is "1"
//  then StepPdf returns zero below x0 and one above (as
//  was true previously); if is value is "0" StepPdf
//  returns one below x0 and zero above.
//  The default of Up is "1" so that this modification
//  should break no existing code that uses it successfully

namespace GooFit {

__device__ auto device_Step(fptype *evt, ParameterContainer &pc) -> fptype {
    int id    = pc.getObservable(0);
    fptype x  = RO_CACHE(evt[id]);
    fptype x0 = pc.getParameter(0);
    int Up    = pc.getConstant(0);
    // it seems that the third argument of incrementIndex
    // is the number of Constants, so change from 0 to 1
    // when adding Up as constant argument in host function
    // may be problem with incrementIndex args, so try "default"
    // mds 211021    pc.incrementIndex(1, 1, 1, 1, 1);
    pc.incrementIndex();

    // note that the value of Up is checked in the host
    // method, so the final "else" should never be
    // executed, but the compiler does not know that
    // and might complain that the method can appear
    // to end without providing a return value
    if(1 == Up) {
        return (x > x0 ? 1 : 0);
    }
    //  this is effectively "else"
    return (x > x0 ? 0 : 1);
}

__device__ device_function_ptr ptr_to_Step = device_Step;
device_function_ptr hptr_to_Step           = device_Step;

__host__ StepPdf::StepPdf(std::string n, Observable _x, Variable x0, int Up)
    : GooPdf("StepPdf", n, _x, x0) {
    registerFunction("ptr_to_Step", ptr_to_Step);
    host_fcn_ptr                       = get_device_symbol_address(ptr_to_Step);
    functionPtrToNameMap[host_fcn_ptr] = "StepPdf";
    registerConstant(Up);

    initialize();
}

__host__ auto StepPdf::integrate(fptype lo, fptype hi) const -> fptype {
    // unsigned int *indices = host_indices + parameters;
    fptype x0 = parametersList[0].getValue();
    int Up    = constantsList[0];
    if(1 == Up) {
        return (hi - x0);
    } else {
        return (x0 - lo);
    }
}

} // namespace GooFit
