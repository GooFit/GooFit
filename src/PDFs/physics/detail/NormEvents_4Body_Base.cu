#include <vector>

#include <mcbooster/Generate.h>
#include <mcbooster/GContainers.h>
#include <mcbooster/Vector4R.h>
#include <mcbooster/EvaluateArray.h>

#include <goofit/PDFs/physics/detail/NormEvents_4Body_Base.h>
#include <goofit/PDFs/physics/detail/Dim5.h>

namespace GooFit {

  __host__ unsigned int NormEvents_4Body_Base::generate4BodyNormEvents(
							      long normSeed,
							      unsigned int numNormEventsToGen,
							      const std::vector<mcbooster::GReal_t>& motherAndDaughterMasses,
							      mcbooster::RealVector_d& norm_M12,
							      mcbooster::RealVector_d& norm_M34,
							      mcbooster::RealVector_d& norm_CosTheta12,
							      mcbooster::RealVector_d& norm_CosTheta34,
							      mcbooster::RealVector_d& norm_phi)
  {
    std::vector<mcbooster::GReal_t> masses(motherAndDaughterMasses.cbegin()+1, motherAndDaughterMasses.cend());
    mcbooster::PhaseSpace phsp(motherAndDaughterMasses[0], masses, numNormEventsToGen, normSeed);
    phsp.Generate(mcbooster::Vector4R(motherAndDaughterMasses[0], 0.0, 0.0, 0.0));
    phsp.Unweight();

    auto nAcc                     = phsp.GetNAccepted();
    mcbooster::BoolVector_d flags = phsp.GetAccRejFlags();
    auto d1                       = phsp.GetDaughters(0);
    auto d2                       = phsp.GetDaughters(1);
    auto d3                       = phsp.GetDaughters(2);
    auto d4                       = phsp.GetDaughters(3);

    auto zip_begin = thrust::make_zip_iterator(thrust::make_tuple(d1.begin(), d2.begin(), d3.begin(), d4.begin()));
    auto zip_end   = zip_begin + d1.size();
    auto new_end   = thrust::remove_if(zip_begin, zip_end, flags.begin(), thrust::logical_not<bool>());

    d1.erase(thrust::get<0>(new_end.get_iterator_tuple()), d1.end());
    d2.erase(thrust::get<1>(new_end.get_iterator_tuple()), d2.end());
    d3.erase(thrust::get<2>(new_end.get_iterator_tuple()), d3.end());
    d4.erase(thrust::get<3>(new_end.get_iterator_tuple()), d4.end());

    mcbooster::ParticlesSet_d pset(4);
    pset[0] = &d1;
    pset[1] = &d2;
    pset[2] = &d3;
    pset[3] = &d4;

    norm_M12.resize(nAcc);
    norm_M34.resize(nAcc);
    norm_CosTheta12.resize(nAcc);
    norm_CosTheta34.resize(nAcc);
    norm_phi.resize(nAcc);
   
    mcbooster::VariableSet_d VarSet(5);
    VarSet[0] = &norm_M12;
    VarSet[1] = &norm_M34;
    VarSet[2] = &norm_CosTheta12;
    VarSet[3] = &norm_CosTheta34;
    VarSet[4] = &norm_phi;

    Dim5 eval = Dim5();
    mcbooster::EvaluateArray<Dim5>(eval, pset, VarSet);

    phsp.FreeResources();
    
    return nAcc;
  }

} // end namespace GooFit
