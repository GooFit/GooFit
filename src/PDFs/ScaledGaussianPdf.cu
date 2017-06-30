#include "goofit/PDFs/basic/ScaledGaussianPdf.h"
#include "goofit/Variable.h"

//#include <limits>

namespace GooFit {

__device__ fptype device_ScaledGaussian(fptype *evt, ParameterContainer &pc) {
    int id = pc.constants[pc.constantIdx + 1];
    fptype x     = evt[id];
    fptype mean  = pc.parameters[pc.parameterIdx + 1] + pc.parameters[pc.parameterIdx + 3];
    fptype sigma = pc.parameters[pc.parameterIdx + 2] * (1 + pc.parameters[pc.parameterIdx + 4]);
    fptype ret   = exp(-0.5 * (x - mean) * (x - mean) / (sigma * sigma));

    pc.incrementIndex (1, 4, 1, 0, 1);

    return ret;
}

__device__ device_function_ptr ptr_to_ScaledGaussian = device_ScaledGaussian;

__host__ ScaledGaussianPdf::ScaledGaussianPdf(
    std::string n, Variable *_x, Variable *mean, Variable *sigma, Variable *delta, Variable *epsilon)
    : GooPdf(_x, n) {

    //placeholder for _x index
    constantsList.push_back (0);

    registerParameter(mean);
    registerParameter(sigma);
    registerParameter(delta);
    registerParameter(epsilon);

    std::vector<unsigned int> pindices;
    //pindices.push_back(mean->getIndex());
    //pindices.push_back(sigma->getIndex());
    //pindices.push_back(delta->getIndex());
    //pindices.push_back(epsilon->getIndex());
    GET_FUNCTION_ADDR(ptr_to_ScaledGaussian);
    initialize(pindices);
}

__host__ void ScaledGaussianPdf::recursiveSetIndices () {
    GET_FUNCTION_ADDR(ptr_to_ScaledGaussian);

    GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName (), "ptr_to_ScaledGaussian");
    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx = num_device_functions++;

    populateArrays ();
}

} // namespace GooFit
