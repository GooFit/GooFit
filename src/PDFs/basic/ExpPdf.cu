#include <goofit/Error.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/basic/ExpPdf.h>

namespace GooFit {

__device__ auto device_Exp(fptype *evt, ParameterContainer &pc) -> fptype {
    int id       = pc.getObservable(0);
    fptype alpha = pc.getParameter(0);
    fptype x     = RO_CACHE(evt[id]);

    fptype ret = exp(alpha * x);

    pc.incrementIndex(1, 1, 0, 1, 1);

    return ret;
}

__device__ auto device_ExpOffset(fptype *evt, ParameterContainer &pc) -> fptype {
    int id   = pc.getObservable(0);
    fptype x = RO_CACHE(evt[id]);
    x -= pc.getParameter(0);
    fptype alpha = pc.getParameter(1);

    fptype ret = exp(alpha * x);

    pc.incrementIndex(1, 2, 0, 1, 1);

    return ret;
}

__device__ auto device_ExpPoly(fptype *evt, ParameterContainer &pc) -> fptype {
    int id   = pc.getObservable(0);
    fptype x = RO_CACHE(evt[id]);

    fptype exparg = 0;

    int np = pc.getNumParameters();
    for(int i = 0; i < np; ++i) {
        exparg += pow(x, i) * pc.getParameter(i);
    }

    fptype ret = exp(exparg);

    pc.incrementIndex(1, np, 0, 1, 1);

    return ret;
}

__device__ auto device_ExpPolyOffset(fptype *evt, ParameterContainer &pc) -> fptype {
    int id   = pc.getObservable(0);
    fptype x = RO_CACHE(evt[id]);
    x -= pc.getParameter(0);

    fptype exparg = 0;

    int np = pc.getNumParameters();
    for(int i = 1; i < np; ++i) {
        exparg += pow(x, i) * pc.getParameter(i);
    }

    fptype ret = exp(exparg);

    pc.incrementIndex(1, np, 0, 1, 1);

    return ret;
}

__device__ device_function_ptr ptr_to_Exp           = device_Exp;
__device__ device_function_ptr ptr_to_ExpPoly       = device_ExpPoly;
__device__ device_function_ptr ptr_to_ExpOffset     = device_ExpOffset;
__device__ device_function_ptr ptr_to_ExpPolyOffset = device_ExpPolyOffset;

__host__ ExpPdf::ExpPdf(std::string n, Observable _x, Variable alpha)
    : GooPdf("ExpPdf", n, _x, alpha) {
    registerFunction("ptr_to_Exp", ptr_to_Exp);
    host_fcn_ptr                       = get_device_symbol_address(ptr_to_Exp);
    functionPtrToNameMap[host_fcn_ptr] = "Exp";

    initialize();
}

__host__ ExpPdf::ExpPdf(std::string n, Observable _x, Variable alpha, Variable offset)
    : GooPdf("ExpPdf", n, _x, offset, alpha) {
    registerFunction("ptr_to_ExpOffset", ptr_to_ExpOffset);
    host_fcn_ptr                       = get_device_symbol_address(ptr_to_ExpOffset);
    functionPtrToNameMap[host_fcn_ptr] = "ExpOffset";
    initialize();
}

__host__ ExpPdf::ExpPdf(std::string n, Observable _x, std::vector<Variable> &weights)
    : GooPdf("ExpPdf", n, _x) {
    if(weights.empty())
        throw GooFit::GeneralError("Weights are empty!");

    for(Variable &w : weights)
        registerParameter(w);

    registerFunction("ptr_to_ExpPoly", ptr_to_ExpPoly);
    host_fcn_ptr                       = get_device_symbol_address(ptr_to_ExpPoly);
    functionPtrToNameMap[host_fcn_ptr] = "ExpPoly";

    initialize();
}

__host__ ExpPdf::ExpPdf(std::string n, Observable _x, std::vector<Variable> &weights, Variable offset)
    : GooPdf("ExpPdf", n, _x, offset) {
    if(weights.empty())
        throw GooFit::GeneralError("Weights are empty!");

    for(Variable &w : weights)
        registerParameter(w);

    registerFunction("ptr_to_ExpPolyOffset", ptr_to_ExpPolyOffset);
    host_fcn_ptr                       = get_device_symbol_address(ptr_to_ExpPolyOffset);
    functionPtrToNameMap[host_fcn_ptr] = "ExpPolyOffset";

    initialize();
}

__host__ auto ExpPdf::integrate(fptype lo, fptype hi) const -> fptype {
    fptype alpha = host_parameters[parametersIdx + 1];

    if(0 == alpha) {
        // This gives a constant 1 all across the range
        return (hi - lo);
    }

    fptype ret = exp(alpha * hi) - exp(alpha * lo);
    ret /= alpha;
    return ret;
}

} // namespace GooFit
