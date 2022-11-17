#include <algorithm>
#include <boost/filesystem.hpp>
#include <cctype>
#include <ctime>
#include <functional>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <thrust/count.h>
#include <vector>
#include <fstream>

#include <mcbooster/Evaluate.h>
#include <mcbooster/EvaluateArray.h>
#include <mcbooster/GContainers.h>
#include <mcbooster/GFunctional.h>
#include <mcbooster/GTypes.h>
#include <mcbooster/Generate.h>
#include <mcbooster/Vector4R.h>

#include "ToyModel.h"

#include <goofit/PDFs/physics/detail/NormEvents_4Body_WeightedDevice.h>
#include <goofit/PDFs/physics/detail/Dim5.h>

using namespace GooFit;

struct genExp {
    fptype gamma;
    unsigned int offset;
    unsigned int seed;

    __host__ __device__ genExp(unsigned int c, fptype d, unsigned int _seed)
        : gamma(d)
        , offset(c)
        , seed(_seed){};

    __host__ __device__ auto operator()(unsigned int x) const -> fptype {
        thrust::random::default_random_engine rand(seed);
        thrust::uniform_real_distribution<fptype> dist(0, 1);

        rand.discard(x + offset);

        return -log(dist(rand)) / gamma;
    }
};

struct evalInvExp {
    fptype gamma;

    __host__ __device__ evalInvExp(fptype d)
        : gamma(d){};

    __host__ __device__ fptype operator()(fptype t) const {
        fptype pdf = gamma * exp(-t * gamma);
        return 1. / pdf;
    }
};

int main(int argc, char **argv) {
    const fptype START_MIXING_ERROR         = 0.00001;
    const fptype MIXING_LOWER_LIMIT         = -0.15;
    const fptype MIXING_UPPER_LIMIT         = 0.15;
    const unsigned int MODEL_MC_EVENTS_NORM = 6e6;
    const unsigned int seed                 = 654321;

    const fptype xMixingStart  = 0.005;
    const fptype yMixingStart  = 0.005;
    std::string outputFilename = "tddp4_fit_results.txt";

    // Build MC events for integration.
    printf("Building the MC sample for the normalization integral.\n");
    std::vector<mcbooster::GReal_t> masses{ToyModel::PI_MASS, ToyModel::PI_MASS, ToyModel::K_MASS, ToyModel::PI_MASS};

    const uint nMC = 6000000;
    mcbooster::PhaseSpace phsp(ToyModel::D0_MASS, masses, nMC, 0);
    phsp.SetSeed(seed);
    phsp.Generate(mcbooster::Vector4R(ToyModel::D0_MASS, 0., 0., 0.));

    auto d1      = phsp.GetDaughters(0);
    auto d2      = phsp.GetDaughters(1);
    auto d3      = phsp.GetDaughters(2);
    auto d4      = phsp.GetDaughters(3);
    auto weights = phsp.GetWeights();

    mcbooster::ParticlesSet_d pset(4);
    pset[0] = &d1;
    pset[1] = &d2;
    pset[2] = &d3;
    pset[3] = &d4;

    auto SigGen_M12_d        = mcbooster::RealVector_d(nMC);
    auto SigGen_M34_d        = mcbooster::RealVector_d(nMC);
    auto SigGen_CosTheta12_d = mcbooster::RealVector_d(nMC);
    auto SigGen_CosTheta34_d = mcbooster::RealVector_d(nMC);
    auto SigGen_phi_d        = mcbooster::RealVector_d(nMC);
    auto SigGen_dtime_d      = mcbooster::RealVector_d(nMC);
    auto SigGen_sigma_d      = mcbooster::RealVector_d(nMC);

    thrust::counting_iterator<unsigned int> index_sequence_begin(0);

    mcbooster::VariableSet_d VarSet_d(5);
    VarSet_d[0] = &SigGen_M12_d;
    VarSet_d[1] = &SigGen_M34_d;
    VarSet_d[2] = &SigGen_CosTheta12_d;
    VarSet_d[3] = &SigGen_CosTheta34_d;
    VarSet_d[4] = &SigGen_phi_d;

    Dim5 eval = Dim5();
    mcbooster::EvaluateArray<Dim5>(eval, pset, VarSet_d);

    // Copy back to the host.
    // TODO: Make this optional.
    mcbooster::RealVector_h SigGen_M12_h        = SigGen_M12_d;
    mcbooster::RealVector_h SigGen_M34_h        = SigGen_M34_d;
    mcbooster::RealVector_h SigGen_CosTheta12_h = SigGen_CosTheta12_d;
    mcbooster::RealVector_h SigGen_CosTheta34_h = SigGen_CosTheta34_d;
    mcbooster::RealVector_h SigGen_phi_h        = SigGen_phi_d;

    // Generate the decay time.
    fptype gamma = 1. / 0.41;
    thrust::transform(index_sequence_begin, index_sequence_begin + nMC, SigGen_dtime_d.begin(), genExp(0, gamma, seed));
    mcbooster::RealVector_h SigGen_dtime_h = SigGen_dtime_d;

    // Generate a decay time weight.
    auto dtime_weights_d = mcbooster::RealVector_d(nMC);
    thrust::transform(SigGen_dtime_d.begin(), SigGen_dtime_d.end(), dtime_weights_d.begin(), evalInvExp(gamma));

    // Get the total weight.
    mcbooster::RealVector_d phsp_weights_d = weights;
    mcbooster::RealVector_d total_weights_d(nMC);
    thrust::transform(dtime_weights_d.begin(),
                      dtime_weights_d.end(),
                      phsp_weights_d.begin(),
                      total_weights_d.begin(),
                      thrust::multiplies<float>());
    mcbooster::RealVector_h SigGen_weights_h = total_weights_d;

    // Get a dummy vector for sigma.
    thrust::fill(SigGen_sigma_d.begin(), SigGen_sigma_d.end(), 0.);
    mcbooster::RealVector_h SigGen_sigma_h = SigGen_sigma_d;

    // Build the norm events.
    NormEvents_4Body_WeightedDevice normEvents(SigGen_M12_d,
                                               SigGen_M34_d,
                                               SigGen_CosTheta12_d,
                                               SigGen_CosTheta34_d,
                                               SigGen_phi_d,
                                               SigGen_dtime_d,
                                               SigGen_sigma_d,
                                               total_weights_d);
    std::vector<NormEvents_4Body_Base *> normEventsVector;
    normEventsVector.push_back(dynamic_cast<NormEvents_4Body_Base *>(&normEvents));

    // Build the model.
    printf("Building the model...\n");
    ToyModel model(xMixingStart, yMixingStart, normEventsVector);
    model.setXMixingRangeForFit(START_MIXING_ERROR, MIXING_LOWER_LIMIT, MIXING_UPPER_LIMIT);
    model.setYMixingRangeForFit(START_MIXING_ERROR, MIXING_LOWER_LIMIT, MIXING_UPPER_LIMIT);

    // Read the input data.
    printf("Reading input data...\n");
    uint fitEvts = 0;
    double m12, m34, c12, c34, phi, dtime;
    std::string inputFile = "tddp4_data.txt";
    std::ifstream file(inputFile);
    std::string line;
    while(std::getline(file, line)) {
        std::stringstream ss(line);
        file >> m12 >> m34 >> c12 >> c34 >> phi >> dtime;
        model.addEventToCurrentDataToFit(m12, m34, c12, c34, phi, dtime, 0.0, fitEvts++);
    }
    file.close();

    // Create the MC file to use for plotting later.
    printf("Creating MC dataset for plotting...\n");
    double tm12, tm34, tc12, tc34, tphi, tdtime, tweight, tsigma, tpdf;
    for(unsigned int i = 0; i < SigGen_weights_h.size(); i++) {
        tm12   = SigGen_M12_h[i];
        tm34   = SigGen_M34_h[i];
        tc12   = SigGen_CosTheta12_h[i];
        tc34   = SigGen_CosTheta34_h[i];
        tphi   = SigGen_phi_h[i];
        tdtime = SigGen_dtime_h[i];
        tsigma = SigGen_sigma_h[i];
        model.addEventToMCToPlot(tm12, tm34, tc12, tc34, tphi, tdtime, 0.0, i);
    }

    // Perform the fit.
    std::vector<std::vector<fptype>> mcValues = model.fitCurrentData(0, outputFilename);

    // Create the output MC file.
    std::ofstream mcfile;
    mcfile.open("tddp4_mc.txt");
    for(unsigned int i = 0; i < SigGen_weights_h.size(); i++) {
        tm12    = SigGen_M12_h[i];
        tm34    = SigGen_M34_h[i];
        tc12    = SigGen_CosTheta12_h[i];
        tc34    = SigGen_CosTheta34_h[i];
        tphi    = SigGen_phi_h[i];
        tdtime  = SigGen_dtime_h[i];
        tweight = SigGen_weights_h[i];
        tpdf    = mcValues[0][i];
        mcfile << tm12 << " " << tm34 << " " << tc12 << " " << tc34 << " " << tphi << " " << tdtime << " " << tweight
               << " " << tpdf << std::endl;
    }
    mcfile.close();

    printf("Done!\n");
    return 0;
}
