#include <goofit/PDFs/physics/detail/Dim5.h>

#include <mcbooster/Evaluate.h>
#include <mcbooster/EvaluateArray.h>
#include <mcbooster/GContainers.h>
#include <mcbooster/GFunctional.h>
#include <mcbooster/GTypes.h>
#include <mcbooster/Vector4R.h>

#include <goofit/Catch.h>
#include <goofit/TestUtils.h>

using namespace GooFit;

TEST_CASE("Test operator()")
{
    const fptype COMPARE_EPS = 1e-5;
    
    const bool PRINT_ME = false;

    // the 4 vecs we want to convert
    // E, px, py, pz [GeV]
    const auto k_4Vec = mcbooster::Vector4R(
        0.65905276036464722, 
        -0.22605460233259722, 
        0.37058687639201848,
        -0.046885439376411875);
    const auto osPi1_4Vec = mcbooster::Vector4R(
        0.35908482669738223, 
        0.075397408921232992, 
        0.24469544143911467, 
        0.20952672690121868);
    const auto osPi2_4Vec = mcbooster::Vector4R(
        0.41772611931236503, 
        0.07358860140319394, 
        -0.24208436188963289,
        -0.30165403210059527);
    const auto ssPi_4Vec = mcbooster::Vector4R(
        0.42897629362560541, 
        0.077068592008170317, 
        -0.37319795594150029, 
        0.13901274457578858);

    // expected vals [GeV] [rad from -pi to pi]
    const fptype expectedM12 = 0.7803089370361741;
    const fptype expectedM34 = 0.5377572566735048;
    const fptype expectedCos12 = -0.22021501289790646;
    const fptype expectedCos34 = -0.02012576600161909;
    const fptype expectedPhi = -2.3558812821;
    if (PRINT_ME)
    {
        GOOFIT_INFO("Values");
    }

    // format input
    mcbooster::Particles_d d1(1);
    d1[0] = k_4Vec;
    mcbooster::Particles_d d2(1);
    d2[0] = osPi1_4Vec;
    mcbooster::Particles_d d3(1);
    d3[0] = osPi2_4Vec;
    mcbooster::Particles_d d4(1);
    d4[0] = ssPi_4Vec;
    if (PRINT_ME)
    {
        for(int i = 0; i < d1.size(); i++)
        {
            std::cout << "D[" << i << "] = " << d1[i] << std::endl;
        }
        GOOFIT_INFO("4 vec values")
    }
    mcbooster::ParticlesSet_d pset(4);
    pset[0] = &d1;
    pset[1] = &d2;
    pset[2] = &d3;
    pset[3] = &d4;
    if (PRINT_ME)
    {
        GOOFIT_INFO("Done formatting input");
    }

    // for output
    auto m12_d        = mcbooster::RealVector_d(1);
    auto m34_d        = mcbooster::RealVector_d(1);
    auto cosTheta12_d = mcbooster::RealVector_d(1);
    auto cosTheta34_d = mcbooster::RealVector_d(1);
    auto phi_d        = mcbooster::RealVector_d(1);
    mcbooster::VariableSet_d VarSet_d(5);
    VarSet_d[0] = &m12_d;
    VarSet_d[1] = &m34_d;
    VarSet_d[2] = &cosTheta12_d;
    VarSet_d[3] = &cosTheta34_d;
    VarSet_d[4] = &phi_d;
    if (PRINT_ME)
    {
        GOOFIT_INFO("Done creating output containers");
    }

    // do conversion
    Dim5 eval = Dim5();
    mcbooster::EvaluateArray<Dim5>(eval, pset, VarSet_d);
    if (PRINT_ME)
    {
        GOOFIT_INFO("Conversion done");
    }

    // copy results to host
    auto m12_h        = mcbooster::RealVector_h(m12_d);
    auto m34_h        = mcbooster::RealVector_h(m34_d);
    auto cosTheta12_h = mcbooster::RealVector_h(cosTheta12_d);
    auto cosTheta34_h = mcbooster::RealVector_h(cosTheta34_d);
    auto phi_h        = mcbooster::RealVector_h(phi_d);
    if (PRINT_ME)
    {
        GOOFIT_INFO("Done copying values to host");
        GOOFIT_INFO("m12 = {}", m12_h[0]);
        GOOFIT_INFO("m34 = {}", m34_h[0]);
        GOOFIT_INFO("c12 = {}", cosTheta12_h[0]);
        GOOFIT_INFO("c34 = {}", cosTheta34_h[0]);
        GOOFIT_INFO("phi = {}", phi_h[0]);
    }

    // check values
    TestUtils::testCompareFloat(m12_h[0], expectedM12, COMPARE_EPS, PRINT_ME, "m12");
    TestUtils::testCompareFloat(m34_h[0], expectedM34, COMPARE_EPS, PRINT_ME, "m34");
    TestUtils::testCompareFloat(cosTheta12_h[0], expectedCos12, COMPARE_EPS, PRINT_ME, "cos12");
    TestUtils::testCompareFloat(cosTheta34_h[0], expectedCos34, COMPARE_EPS, PRINT_ME, "cos34");
    TestUtils::testCompareFloat(phi_h[0], expectedPhi, COMPARE_EPS, PRINT_ME, "phi");
}