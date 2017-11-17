#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <goofit/PDFs/physics/EvalVar.h>
#include <goofit/PDFs/physics/UserUtils.h>
#include <mcbooster/EvaluateArray.h>

using namespace GooFit;

TEST(Vectors, Convert5param4bodySpace) {
    Eigen::Matrix<fptype, 10, 16> matI;

    matI << -0.0200374, 0.246687, 0.410008, 0.687807, 0.203659, -0.0540155, 0.125961, 0.282384, -0.30315, -0.14159,
        -0.593475, 0.695442, 0.119529, -0.0510815, 0.057506, 0.199206, -0.533425, 0.24162, -0.34825, 0.841379, 0.467861,
        -0.190084, 0.0466843, 0.526008, 0.114051, 0.00314688, 0.287416, 0.339272, -0.0484866, -0.0546832, 0.0141501,
        0.158181, -0.0825251, 0.103312, 0.273959, 0.579875, 0.156211, 0.386267, -0.193338, 0.480066, 0.162955, -0.25586,
        0.16105, 0.370723, -0.236641, -0.233719, -0.241672, 0.434176, 0.414422, 0.0156387, 0.118724, 0.655594,
        0.0188908, -0.198626, 0.346841, 0.423778, -0.0527912, -0.0564821, -0.0895424, 0.182961, -0.380522, 0.239469,
        -0.376023, 0.602508, -0.052077, -0.608464, 0.127396, 0.795542, 0.208922, 0.317326, 0.220122, 0.460736,
        0.0520564, 0.0696663, -0.122869, 0.20528, -0.208901, 0.221471, -0.224649, 0.403282, -0.185982, -0.218172,
        -0.132677, 0.586096, -0.100874, 0.0636139, -0.0320129, 0.186351, -0.365637, -0.00680145, 0.0279881, 0.392428,
        0.652492, 0.16136, 0.136702, 0.699965, -0.190126, -0.114267, -0.0734665, 0.546186, -0.173961, 0.153112,
        0.308143, 0.410046, 0.553302, 0.133658, -0.162681, 0.608237, -0.189214, -0.172503, -0.0719956, 0.300371,
        -0.476565, -0.376938, 0.121717, 0.792294, -0.0319598, 0.0709594, 0.0715396, 0.175084, 0.297197, 0.271234,
        0.165535, 0.456919, 0.211328, 0.0347447, -0.358791, 0.440543, 0.00358106, -0.05433, -0.0144613, 0.496881,
        -0.323238, -0.0803101, 0.0576479, 0.365699, -0.273376, 0.0620731, 0.165497, 0.354199, 0.593034, 0.072567,
        -0.208683, 0.648061, 0.118637, 0.115931, 0.270851, 0.58702, 0.431136, -0.0040225, -0.314793, 0.551787,
        -0.0572405, 0.0658558, -0.0581859, 0.174582, -0.492533, -0.177765, 0.102128, 0.551451;

    Eigen::Matrix<fptype, 10, 5> finalsI;

    finalsI << 0.763655, 0.665051, 0.109243, 0.913621, -3.07925, 1.33111, 0.386733, 0.536461, 0.836249, 1.16468,
        0.933733, 0.62948, -0.506615, -0.141132, 1.98515, 0.852672, 0.423044, -0.250731, -0.843664, -1.35136, 1.16108,
        0.374457, 0.328627, -0.61922, 1.3365, 0.680718, 1.02958, 0.226766, -0.875116, 0.818193, 0.851609, 0.797763,
        -0.60318, 0.755663, -0.194568, 0.739088, 0.644879, 0.779099, 0.0291094, 0.552592, 0.788588, 0.939335, -0.96629,
        -0.880504, 0.383515, 0.990039, 0.458714, -0.557693, -0.843915, -1.20109;

    Eigen::Matrix<fptype, 16, 10> mat   = matI.transpose();
    Eigen::Matrix<fptype, 5, 10> finals = finalsI.transpose();

    mcbooster::Particles_h h_p0(10);
    mcbooster::Particles_h h_p1(10);
    mcbooster::Particles_h h_p2(10);
    mcbooster::Particles_h h_p3(10);

    for(int i = 0; i < mat.cols(); i++) {
        h_p0[i].set(mat(3, i), mat(0, i), mat(1, i), mat(2, i));
        h_p1[i].set(mat(7, i), mat(4, i), mat(5, i), mat(6, i));
        h_p2[i].set(mat(11, i), mat(8, i), mat(9, i), mat(10, i));
        h_p3[i].set(mat(15, i), mat(12, i), mat(13, i), mat(14, i));
    }

    mcbooster::Particles_d d_p0(h_p0.begin(), h_p0.end());
    mcbooster::Particles_d d_p1(h_p1.begin(), h_p1.end());
    mcbooster::Particles_d d_p2(h_p2.begin(), h_p2.end());
    mcbooster::Particles_d d_p3(h_p3.begin(), h_p3.end());

    mcbooster::RealVector_d d_m12(10);
    mcbooster::RealVector_d d_m34(10);
    mcbooster::RealVector_d d_cos12(10);
    mcbooster::RealVector_d d_cos34(10);
    mcbooster::RealVector_d d_phi(10);

    mcbooster::ParticlesSet_d particles = {&d_p0, &d_p1, &d_p2, &d_p3};
    mcbooster::VariableSet_d variables  = {&d_m12, &d_m34, &d_cos12, &d_cos34, &d_phi};

    Dim5 eval = Dim5();
    mcbooster::EvaluateArray<Dim5>(eval, particles, variables);

    mcbooster::RealVector_h h_m12(d_m12.begin(), d_m12.end());
    mcbooster::RealVector_h h_m34(d_m34.begin(), d_m34.end());
    mcbooster::RealVector_h h_cos12(d_cos12.begin(), d_cos12.end());
    mcbooster::RealVector_h h_cos34(d_cos34.begin(), d_cos34.end());
    mcbooster::RealVector_h h_phi(d_phi.begin(), d_phi.end());

    for(int i = 0; i < mat.cols(); i++) {
        SCOPED_TRACE("i = " + std::to_string(i));
        ASSERT_NEAR(h_m12[i], finals(0, i), .00001);
        ASSERT_NEAR(h_m34[i], finals(1, i), .00001);
        ASSERT_NEAR(h_cos12[i], finals(2, i), .00001);
        ASSERT_NEAR(h_cos34[i], finals(3, i), .00001);
        ASSERT_NEAR(h_phi[i], finals(4, i), .00001);
    }

    auto new_finals = to_5param(mat);
    for(int i = 0; i < 5; i++) {
        for(int j = 0; j < mat.cols(); j++) {
            SCOPED_TRACE("i = " + std::to_string(i) + ", j = " + std::to_string(j));
            ASSERT_NEAR(finals(i, j), new_finals(i, j), .00001);
        }
    }
}
