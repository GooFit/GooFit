#include "goofit/PDFs/physics/DalitzVetoPdf.h"
#include "goofit/PDFs/physics/DalitzPlotHelpers.h"

namespace GooFit {

__device__ fptype device_DalitzVeto(fptype *evt, ParameterContainer &pc) {
    int numConstants   = RO_CACHE(pc.constants[pc.constantIdx]);
    int numObservables = RO_CACHE(pc.observables[pc.observableIdx]);

    int idx1 = RO_CACHE(pc.observables[pc.observableIdx + 1]);
    int idx2 = RO_CACHE(pc.observables[pc.observableIdx + 2]);

    fptype motherM = RO_CACHE(pc.parameters[pc.parameterIdx + 1]);
    fptype d1m     = RO_CACHE(pc.parameters[pc.parameterIdx + 2]);
    fptype d2m     = RO_CACHE(pc.parameters[pc.parameterIdx + 3]);
    fptype d3m     = RO_CACHE(pc.parameters[pc.parameterIdx + 4]);

    fptype x = evt[idx1];
    fptype y = evt[idx2];

    fptype massSum = motherM * motherM + d1m * d1m + d2m * d2m + d3m * d3m;
    fptype z       = massSum - x - y;

    fptype ret            = inDalitz(x, y, motherM, d1m, d2m, d3m) ? 1.0 : 0.0;
    unsigned int numVetos = RO_CACHE(pc.constants[pc.constantIdx + 1]);

    for(int i = 0; i < numVetos; ++i) {
        unsigned int varIndex = pc.constants[pc.constantIdx + 2 + i];
        fptype minimum        = RO_CACHE(pc.parameters[pc.parameterIdx + 5 + i * 2]);
        fptype maximum        = RO_CACHE(pc.parameters[pc.parameterIdx + 5 + i * 2 + 1]);
        fptype currDalitzVar  = (PAIR_12 == varIndex ? x : PAIR_13 == varIndex ? y : z);

        ret *= ((currDalitzVar < maximum) && (currDalitzVar > minimum)) ? 0.0 : 1.0;
    }

    // TODO: Prefer this function, not incrementIndex();
    // pc.incrementIndex(1, numVetos*2 + 4, numConstants, numObservables, 1);
    pc.incrementIndex();
    return ret;
}

__device__ device_function_ptr ptr_to_DalitzVeto = device_DalitzVeto;

__host__ DalitzVetoPdf::DalitzVetoPdf(std::string n,
                                      Observable _x,
                                      Observable _y,
                                      Variable motherM,
                                      Variable d1m,
                                      Variable d2m,
                                      Variable d3m,
                                      std::vector<VetoInfo> vetos)
    : GooPdf(n, _x, _y) {
    std::vector<unsigned int> pindices;
    pindices.push_back(registerParameter(motherM));
    pindices.push_back(registerParameter(d1m));
    pindices.push_back(registerParameter(d2m));
    pindices.push_back(registerParameter(d3m));

    pindices.push_back(vetos.size());
    constantsList.push_back(vetos.size());

    for(auto &veto : vetos) {
        pindices.push_back(veto.cyclic_index);
        pindices.push_back(registerParameter(veto.minimum));
        pindices.push_back(registerParameter(veto.maximum));
        constantsList.push_back(veto.cyclic_index);
    }

    GET_FUNCTION_ADDR(ptr_to_DalitzVeto);
    initialize(pindices);
}

void DalitzVetoPdf::recursiveSetIndices() {
    GET_FUNCTION_ADDR(ptr_to_DalitzVeto);

    GOOFIT_TRACE("host_function_table[{}]({})", num_device_functions, getName(), "ptr_to_DalitzVeto");
    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx                               = num_device_functions++;

    populateArrays();
}
} // namespace GooFit
