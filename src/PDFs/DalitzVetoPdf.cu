#include "goofit/PDFs/physics/DalitzVetoPdf.h"
#include "goofit/PDFs/physics/DalitzPlotHelpers.h"

namespace GooFit {

__device__ fptype device_DalitzVeto(fptype *evt, ParameterContainer &pc) {
    int numConstants = RO_CACHE(pc.constants[pc.constantIdx]);
    int idx1 = RO_CACHE(pc.constants[pc.constantIdx + 2]);
    int idx2 = RO_CACHE(pc.constants[pc.constantIdx + 3]);

    fptype motherM = RO_CACHE(pc.parameters[pc.parameterIdx + 1]);
    fptype d1m     = RO_CACHE(pc.parameters[pc.parameterIdx + 2]);
    fptype d2m     = RO_CACHE(pc.parameters[pc.parameterIdx + 3]);
    fptype d3m     = RO_CACHE(pc.parameters[pc.parameterIdx + 4]);

    fptype x = evt[idx1];
    fptype y = evt[idx2];

    fptype massSum = motherM * motherM + d1m * d1m + d2m * d2m + d3m * d3m;
    fptype z       = massSum - x - y;

    fptype ret            = inDalitz(x, y, motherM, d1m, d2m, d3m) ? 1.0 : 0.0;
    unsigned int numVetos = RO_CACHE(pc.constants[pc.constantIdx + 4]);

    for(int i = 0; i < numVetos; ++i) {
        unsigned int varIndex = pc.constants[pc.constantIdx + 5 + i];
        fptype minimum        = RO_CACHE(pc.parameters[pc.parameterIdx + 5 + i * 2]);
        fptype maximum        = RO_CACHE(pc.parameters[pc.parameterIdx + 5 + i * 2 + 1]);
        fptype currDalitzVar  = (PAIR_12 == varIndex ? x : PAIR_13 == varIndex ? y : z);

        ret *= ((currDalitzVar < maximum) && (currDalitzVar > minimum)) ? 0.0 : 1.0;
    }

    pc.incrementIndex(1, 4, numConstants, 0, 1);

    return ret;
}

__device__ device_function_ptr ptr_to_DalitzVeto = device_DalitzVeto;

__host__ DalitzVetoPdf::DalitzVetoPdf(std::string n,
                                      Variable *_x,
                                      Variable *_y,
                                      Variable *motherM,
                                      Variable *d1m,
                                      Variable *d2m,
                                      Variable *d3m,
                                      std::vector<VetoInfo *> vetos)
    : GooPdf(nullptr, n) {

    registerObservable(_x);
    registerObservable(_y);

    constantsList.push_back(observablesList.size());
    constantsList.push_back(0);
    constantsList.push_back(0);

    std::vector<unsigned int> pindices;
    pindices.push_back(registerParameter(motherM));
    pindices.push_back(registerParameter(d1m));
    pindices.push_back(registerParameter(d2m));
    pindices.push_back(registerParameter(d3m));

    pindices.push_back(vetos.size());
    constantsList.push_back (vetos.size());

    for(auto &veto : vetos) {
        pindices.push_back(veto->cyclic_index);
        pindices.push_back(registerParameter(veto->minimum));
        pindices.push_back(registerParameter(veto->maximum));
        constantsList.push_back (veto->cyclic_index);
    }

    GET_FUNCTION_ADDR(ptr_to_DalitzVeto);
    initialize(pindices);
}

void DalitzVetoPdf::recursiveSetIndices () {
}

} // namespace GooFit
