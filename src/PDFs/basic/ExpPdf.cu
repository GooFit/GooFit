#include <goofit/Error.h>
#include <goofit/PDFs/basic/ExpPdf.h>

namespace GooFit {

__device__ fptype device_Exp(fptype *evt, fptype *p, unsigned int *indices) {
    fptype x     = evt[indices[2 + indices[0]]];
    fptype alpha = p[indices[1]];

    fptype ret = exp(alpha * x);
    return ret;
}

__device__ fptype device_ExpOffset(fptype *evt, fptype *p, unsigned int *indices) {
    fptype x = evt[indices[2 + indices[0]]];
    x -= p[indices[1]];
    fptype alpha = p[indices[2]];

    fptype ret = exp(alpha * x);
    return ret;
}

__device__ fptype device_ExpPoly(fptype *evt, fptype *p, unsigned int *indices) {
    fptype x = evt[indices[2 + indices[0]]];

    fptype exparg = 0;

    for(int i = 0; i <= indices[0]; ++i) {
        exparg += pow(x, i) * p[indices[i + 1]];
    }

    fptype ret = exp(exparg);
    return ret;
}

__device__ fptype device_ExpPolyOffset(fptype *evt, fptype *p, unsigned int *indices) {
    fptype x = evt[indices[2 + indices[0]]];
    x -= p[indices[1]];

    fptype exparg = 0;

    for(int i = 0; i <= indices[0]; ++i) {
        exparg += pow(x, i) * p[indices[i + 2]];
    }

    fptype ret = exp(exparg);
    return ret;
}

__device__ device_function_ptr ptr_to_Exp           = device_Exp;
__device__ device_function_ptr ptr_to_ExpPoly       = device_ExpPoly;
__device__ device_function_ptr ptr_to_ExpOffset     = device_ExpOffset;
__device__ device_function_ptr ptr_to_ExpPolyOffset = device_ExpPolyOffset;

__host__ ExpPdf::ExpPdf(std::string n, Observable _x, Variable alpha)
    : GooPdf(n, _x) {
    std::vector<unsigned int> pindices;

    pindices.push_back(registerParameter(alpha));
    GET_FUNCTION_ADDR(ptr_to_Exp);
    initialize(pindices);
}

__host__ ExpPdf::ExpPdf(std::string n, Observable _x, Variable alpha, Variable offset)
    : GooPdf(n, _x) {
    std::vector<unsigned int> pindices;

    pindices.push_back(registerParameter(offset));
    pindices.push_back(registerParameter(alpha));
    GET_FUNCTION_ADDR(ptr_to_ExpOffset);
    initialize(pindices);
}

__host__ ExpPdf::ExpPdf(std::string n, Observable _x, std::vector<Variable> &weights)
    : GooPdf(n, _x) {
    std::vector<unsigned int> pindices;

    if(weights.empty())
        throw GooFit::GeneralError("Weights are empty!");

    for(Variable &w : weights)
        pindices.push_back(registerParameter(w));

    GET_FUNCTION_ADDR(ptr_to_ExpPoly);

    initialize(pindices);
}

__host__ ExpPdf::ExpPdf(std::string n, Observable _x, std::vector<Variable> &weights, Variable offset)
    : GooPdf(n, _x) {
    std::vector<unsigned int> pindices;

    pindices.push_back(registerParameter(offset));

    if(weights.empty())
        throw GooFit::GeneralError("Weights are empty!");

    for(Variable &w : weights)
        pindices.push_back(registerParameter(w));

    GET_FUNCTION_ADDR(ptr_to_ExpPolyOffset);

    initialize(pindices);
}

__host__ fptype ExpPdf::integrate(fptype lo, fptype hi) const {
    fptype alpha = host_params[host_indices[parameters + 1]];

    if(0 == alpha) {
        // This gives a constant 1 all across the range
        return (hi - lo);
    }

    fptype ret = exp(alpha * hi) - exp(alpha * lo);
    ret /= alpha;
    return ret;
}

} // namespace GooFit
