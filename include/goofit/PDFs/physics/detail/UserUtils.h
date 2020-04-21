#pragma once

#include <goofit/PDFs/physics/detail/Dim5.h>

#include <Eigen/Dense>
#include <mcbooster/EvaluateArray.h>

namespace GooFit {

/// This is a helper to convert values. It copies back and forth, so is only to be used in scripts
Eigen::Matrix<fptype, 5, Eigen::Dynamic> to_5param(const Eigen::Matrix<fptype, 16, Eigen::Dynamic> &mat) {
    mcbooster::Particles_h h_p0(mat.cols());
    mcbooster::Particles_h h_p1(mat.cols());
    mcbooster::Particles_h h_p2(mat.cols());
    mcbooster::Particles_h h_p3(mat.cols());

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

    mcbooster::RealVector_d d_m12(mat.cols());
    mcbooster::RealVector_d d_m34(mat.cols());
    mcbooster::RealVector_d d_cos12(mat.cols());
    mcbooster::RealVector_d d_cos34(mat.cols());
    mcbooster::RealVector_d d_phi(mat.cols());

    mcbooster::ParticlesSet_d particles = {&d_p0, &d_p1, &d_p2, &d_p3};
    mcbooster::VariableSet_d variables  = {&d_m12, &d_m34, &d_cos12, &d_cos34, &d_phi};

    Dim5 eval = Dim5();
    mcbooster::EvaluateArray<Dim5>(eval, particles, variables);

    mcbooster::RealVector_h h_m12(d_m12.begin(), d_m12.end());
    mcbooster::RealVector_h h_m34(d_m34.begin(), d_m34.end());
    mcbooster::RealVector_h h_cos12(d_cos12.begin(), d_cos12.end());
    mcbooster::RealVector_h h_cos34(d_cos34.begin(), d_cos34.end());
    mcbooster::RealVector_h h_phi(d_phi.begin(), d_phi.end());

    Eigen::Matrix<fptype, 5, Eigen::Dynamic> output(5, mat.cols());

    for(int i = 0; i < mat.cols(); i++) {
        output(0, i) = h_m12[i];
        output(1, i) = h_m34[i];
        output(2, i) = h_cos12[i];
        output(3, i) = h_cos34[i];
        output(4, i) = h_phi[i];
    }

    return output;
}

} // namespace GooFit
