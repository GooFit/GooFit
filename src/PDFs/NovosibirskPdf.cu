#include "goofit/PDFs/basic/NovosibirskPdf.h"

namespace GooFit {

__device__ fptype device_Novosibirsk(fptype* evt, ParameterContainer &pc) {
    fptype _Mean  = pc.parameters[1];
    fptype _Sigma = pc.parameters[2];
    fptype _Tail  = pc.parameters[3];
    fptype x      = evt[0];

    pc.incrementIndex (1, 3, 0, 0, 1);

    fptype qa = 0;
    fptype qb = 0;
    fptype qc = 0;
    fptype qx = 0;
    fptype qy = 0;

    if(fabs(_Tail) < 1.e-7) {
        qc = 0.5 * POW2((x - _Mean) / _Sigma);
    } else {
        qa = _Tail * sqrt(log(4.));
        qb = sinh(qa) / qa;
        qx = (x - _Mean) / _Sigma * qb;
        qy = 1. + _Tail * qx;

        //---- Cutting curve from right side

        if(qy > 1.e-7)
            qc = 0.5 * (POW2(log(qy) / _Tail) + _Tail * _Tail);
        else
            qc = 15.0;
    }

    //---- Normalize the result
    return exp(-qc);
}

__device__ device_function_ptr ptr_to_Novosibirsk = device_Novosibirsk;

__host__ NovosibirskPdf::NovosibirskPdf(std::string n, Variable *_x, Variable *mean, Variable *sigma, Variable *tail)
    : GooPdf(_x, n) {
    std::vector<unsigned int> pindices;
    pindices.push_back(registerParameter(mean));
    pindices.push_back(registerParameter(sigma));
    pindices.push_back(registerParameter(tail));
    GET_FUNCTION_ADDR(ptr_to_Novosibirsk);
    initialize(pindices);
}

__host__ void NovosibirskPdf::recursiveSetIndices () {
    GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName (), "ptr_to_Novosibirsk");
    GET_FUNCTION_ADDR(ptr_to_Novosibirsk);

    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx = num_device_functions ++;

    populateArrays ();
}

} // namespace GooFit
