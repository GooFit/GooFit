#include <stdio.h>

#include <goofit/PDFs/physics/detail/SFCalculator_TD.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/PDFs/physics/detail/Dim5.h>

namespace GooFit {

SFCalculator_TD::SFCalculator_TD() = default;

__device__ auto SFCalculator_TD::operator()(thrust::tuple<int, fptype *, int> t) const -> fpcomplex {
    // int tid = blockIdx.x*blockDim.x + threadIdx.x;
    // bool printStatus = tid == 0;
    const bool PRINT_ME = false;

    int evtNum  = thrust::get<0>(t);
    
    fptype *evt = thrust::get<1>(t) + (evtNum * thrust::get<2>(t));

    // unsigned int *indices = paramIndices + _parameters; // Jump to DALITZPLOT position within parameters array
    // int parameter_i       = 12 + (2 * indices[6]) + (indices[3] * 2)
    //                  + (_spinfactor_i * 2); // Find position of this resonance relative to DALITZPLOT start
    // unsigned int functn_i = indices[parameter_i];
    // unsigned int params_i = indices[parameter_i + 1];

    ParameterContainer pc;

    // Increment to TDDP function
    while(pc.funcIdx < dalitzFuncId)
        pc.incrementIndex();

    if (PRINT_ME)
    {
        printf("dalitzFuncId: %u, _spinfactor_i: %u, event: %d\n", dalitzFuncId, _spinfactor_i, evtNum);
    }

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

    if (PRINT_ME)
    {
        printf("m12 = %f\n", m12);
        printf("m34 = %f\n", m34);
        printf("cos12 = %f\n", cos12);
        printf("cos34 = %f\n", cos34);
        printf("phi = %f\n", phi);
    }

    fptype M  = pc.getConstant(0);
    fptype m1 = pc.getConstant(1);
    fptype m2 = pc.getConstant(2);
    fptype m3 = pc.getConstant(3);
    fptype m4 = pc.getConstant(4);
    if (PRINT_ME)
    {
        printf("M = %f\n", M);
        printf("m1 = %f\n", m1);
        printf("m2 = %f\n", m2);
        printf("m3 = %f\n", m3);
        printf("m4 = %f\n", m4);
    }

    fptype vecs[16];
    get4Vecs(vecs, m12, m34, cos12, cos34, phi, M, m1, m2, m3, m4);
    // printf("%i, %i, %f, %f, %f, %f, %f \n",evtNum, thrust::get<2>(t), m12, m34, cos12, cos34, phi );
    // printf("vec%i %f, %f, %f, %f\n",0, vecs[0], vecs[1], vecs[2], vecs[3]);
    // printf("vec%i %f, %f, %f, %f\n",1, vecs[4], vecs[5], vecs[6], vecs[7]);
    // printf("vec%i %f, %f, %f, %f\n",2, vecs[8], vecs[9], vecs[10], vecs[11]);
    // printf("vec%i %f, %f, %f, %f\n",3, vecs[12], vecs[13], vecs[14], vecs[15]);

    // increment to spin factor:
    while(pc.funcIdx < _spinfactor_i)
        pc.incrementIndex();

    auto func = reinterpret_cast<spin_function_ptr>(d_function_table[pc.funcIdx]);
    fptype sf = (*func)(vecs, pc);
    // printf("SpinFactors %i : %.7g\n",_spinfactor_i, sf );
    if (PRINT_ME)
    {
        printf("SF value: %f\n", sf);
    }
    return {sf, 0.0};
}

} // namespace GooFit
