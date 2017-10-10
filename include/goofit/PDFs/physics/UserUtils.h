#pragma once

#include "EvalVar.h"

#include <Eigen/Dense>
#include <mcbooster/EvaluateArray.h>

namespace GooFit {
    
/// This is a helper to convert values. It copies back and forth, so is only to be used in scripts
Eigen::Matrix<fptype, Eigen::Dynamic, 5> to_5param(const Eigen::Matrix<fptype, Eigen::Dynamic, 16> &mat) {
    
    mcbooster::Particles_h h_p0(mat.rows());
    mcbooster::Particles_h h_p1(mat.rows());
    mcbooster::Particles_h h_p2(mat.rows());
    mcbooster::Particles_h h_p3(mat.rows());
    
    for(int i=0; i<mat.rows(); i++) {
        h_p0[i].set(mat(i, 3), mat(i, 0), mat(i, 1), mat(i, 2));
        h_p1[i].set(mat(i, 7), mat(i, 4), mat(i, 5), mat(i, 6));
        h_p2[i].set(mat(i,11), mat(i, 8), mat(i, 9), mat(i,10));
        h_p3[i].set(mat(i,15), mat(i,12), mat(i,13), mat(i,14));
    }
    
    mcbooster::Particles_d d_p0(h_p0.begin(), h_p0.end());
    mcbooster::Particles_d d_p1(h_p1.begin(), h_p1.end());
    mcbooster::Particles_d d_p2(h_p2.begin(), h_p2.end());
    mcbooster::Particles_d d_p3(h_p3.begin(), h_p3.end());
    
    mcbooster::RealVector_d d_m12(mat.rows());
    mcbooster::RealVector_d d_m34(mat.rows());
    mcbooster::RealVector_d d_cos12(mat.rows());
    mcbooster::RealVector_d d_cos34(mat.rows());
    mcbooster::RealVector_d d_phi(mat.rows());
    
    mcbooster::ParticlesSet_d particles = {&d_p0, &d_p1, &d_p2, &d_p3};
    mcbooster::VariableSet_d variables = {&d_m12, &d_m34, &d_cos12, &d_cos34, &d_phi};
    
    Dim5 eval = Dim5();
    mcbooster::EvaluateArray<Dim5>(eval, particles, variables);
    
    mcbooster::RealVector_h h_m12(d_m12.begin(), d_m12.end());
    mcbooster::RealVector_h h_m34(d_m34.begin(), d_m34.end());
    mcbooster::RealVector_h h_cos12(d_cos12.begin(), d_cos12.end());
    mcbooster::RealVector_h h_cos34(d_cos34.begin(), d_cos34.end());
    mcbooster::RealVector_h h_phi(d_phi.begin(), d_phi.end());
    
    Eigen::Matrix<fptype, Eigen::Dynamic, 5> output(mat.rows(), 5);
    
    for(int i=0; i<10; i++) {
        output(i,0) = h_m12[i];
        output(i,1) = h_m34[i];
        output(i,2) = h_cos12[i];
        output(i,3) = h_cos34[i];
        output(i,4) = h_phi[i];
    }
    
    return output;
}

} // namespace GooFit
