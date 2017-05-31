#include "goofit/PDFs/physics/DalitzVetoPdf.h"
#include "goofit/PDFs/physics/DalitzPlotHelpers.h"

namespace GooFit {

__device__ fptype device_DalitzVeto(fptype *evt, fptype *p, unsigned int *indices) {
    fptype x = evt[RO_CACHE(indices[2 + RO_CACHE(indices[0]) + 0])];
    fptype y = evt[RO_CACHE(indices[2 + RO_CACHE(indices[0]) + 1])];

    fptype motherM = RO_CACHE(p[RO_CACHE(indices[1])]);
    fptype d1m     = RO_CACHE(p[RO_CACHE(indices[2])]);
    fptype d2m     = RO_CACHE(p[RO_CACHE(indices[3])]);
    fptype d3m     = RO_CACHE(p[RO_CACHE(indices[4])]);

    fptype massSum = motherM * motherM + d1m * d1m + d2m * d2m + d3m * d3m;
    fptype z       = massSum - x - y;

    fptype ret            = inDalitz(x, y, motherM, d1m, d2m, d3m) ? 1.0 : 0.0;
    unsigned int numVetos = RO_CACHE(indices[5]);

    for(int i = 0; i < numVetos; ++i) {
        unsigned int varIndex = indices[6 + i * 3 + 0];
        fptype minimum        = RO_CACHE(p[RO_CACHE(indices[6 + i * 3 + 1])]);
        fptype maximum        = RO_CACHE(p[RO_CACHE(indices[6 + i * 3 + 2])]);
        fptype currDalitzVar  = (PAIR_12 == varIndex ? x : PAIR_13 == varIndex ? y : z);

        ret *= ((currDalitzVar < maximum) && (currDalitzVar > minimum)) ? 0.0 : 1.0;
    }

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

    std::vector<unsigned int> pindices;
    pindices.push_back(registerParameter(motherM));
    pindices.push_back(registerParameter(d1m));
    pindices.push_back(registerParameter(d2m));
    pindices.push_back(registerParameter(d3m));

    pindices.push_back(vetos.size());

    for(auto &veto : vetos) {
        pindices.push_back(veto->cyclic_index);
        pindices.push_back(registerParameter(veto->minimum));
        pindices.push_back(registerParameter(veto->maximum));
    }

    GET_FUNCTION_ADDR(ptr_to_DalitzVeto);
    initialize(pindices);
}
} // namespace GooFit
