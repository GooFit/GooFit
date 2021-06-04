#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/DalitzVetoPdf.h>

namespace GooFit {

__device__ auto device_DalitzVeto(fptype *evt, ParameterContainer &pc) -> fptype {
    int idx1 = pc.getObservable(0);
    int idx2 = pc.getObservable(1);

    fptype motherM = pc.getParameter(0);
    fptype d1m     = pc.getParameter(1);
    fptype d2m     = pc.getParameter(2);
    fptype d3m     = pc.getParameter(3);

    fptype x = RO_CACHE(evt[idx1]);
    fptype y = RO_CACHE(evt[idx2]);

    fptype massSum = motherM * motherM + d1m * d1m + d2m * d2m + d3m * d3m;
    fptype z       = massSum - x - y;

    fptype ret            = inDalitz(x, y, motherM, d1m, d2m, d3m) ? 1.0 : 0.0;
    unsigned int numVetos = pc.getConstant(0);

    for(int i = 0; i < numVetos; ++i) {
        unsigned int varIndex = pc.getConstant(1 + i);
        fptype minimum        = pc.getParameter(4 + i * 2);
        fptype maximum        = pc.getParameter(4 + i * 2 + 1);
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
    : GooPdf("DalitzVetoPdf", n, _x, _y, motherM, d1m, d2m, d3m) {
    registerConstant(vetos.size());

    for(auto &veto : vetos) {
        registerParameter(veto.minimum);
        registerParameter(veto.maximum);
        registerConstant(veto.cyclic_index);
    }

    registerFunction("ptr_to_DalitzVeto", ptr_to_DalitzVeto);

    initialize();
}
} // namespace GooFit
