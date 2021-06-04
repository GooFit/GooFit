#include <goofit/PDFs/physics/detail/SFCalculator.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/PDFs/physics/detail/Dim5.h>

namespace GooFit {

SFCalculator::SFCalculator() = default;

__device__ auto SFCalculator::operator()(thrust::tuple<int, fptype *, int> t) const -> fpcomplex {
    int evtNum  = thrust::get<0>(t);
    fptype *evt = thrust::get<1>(t) + (evtNum * thrust::get<2>(t));

    ParameterContainer pc;

    // Increment to DP
    while(pc.funcIdx < dalitzFuncId)
        pc.incrementIndex();

    int id_m12   = pc.getObservable(0);
    int id_m34   = pc.getObservable(1);
    int id_cos12 = pc.getObservable(2);
    int id_cos34 = pc.getObservable(3);
    int id_phi   = pc.getObservable(4);

    fptype m12   = evt[id_m12];
    fptype m34   = evt[id_m34];
    fptype cos12 = evt[id_cos12];
    fptype cos34 = evt[id_cos34];
    fptype phi   = evt[id_phi];

    fptype M  = pc.getConstant(0);
    fptype m1 = pc.getConstant(1);
    fptype m2 = pc.getConstant(2);
    fptype m3 = pc.getConstant(3);
    fptype m4 = pc.getConstant(4);

    fptype vecs[16];
    get4Vecs(vecs, m12, m34, cos12, cos34, phi, M, m1, m2, m3, m4);
    // printf("%i, %i, %f, %f, %f, %f, %f \n",evtNum, thrust::get<2>(t), m12, m34, cos12, cos34, phi );
    // printf("vec%i %f, %f, %f, %f\n",0, vecs[0], vecs[1], vecs[2], vecs[3]);
    // printf("vec%i %f, %f, %f, %f\n",1, vecs[4], vecs[5], vecs[6], vecs[7]);
    // printf("vec%i %f, %f, %f, %f\n",2, vecs[8], vecs[9], vecs[10], vecs[11]);
    // printf("vec%i %f, %f, %f, %f\n",3, vecs[12], vecs[13], vecs[14], vecs[15]);

    // loop until our appropriate spin factor
    while(pc.funcIdx < _spinfactor_i)
        pc.incrementIndex();

    auto func = reinterpret_cast<spin_function_ptr>(d_function_table[pc.funcIdx]);
    fptype sf = (*func)(vecs, pc);
    // printf("SpinFactors %i : %.7g\n",evtNum, sf );
    return {sf, 0.0};
}

} // namespace GooFit
