// ROOT
#include <TFile.h>
#include <TTree.h>

// GooFit stuff
#include <fstream>
#include <goofit/Application.h>
#include <goofit/PDFs/basic/PolynomialPdf.h>
#include <goofit/PDFs/combine/AddPdf.h>
#include <goofit/PDFs/physics/DP4Pdf.h>
#include <goofit/PDFs/physics/Tddp4Pdf.h>
#include <goofit/PDFs/physics/TruthResolution_Aux.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>
#include <thrust/count.h>

using namespace std;
using namespace GooFit;

// Constants used in more than one PDF component.
const fptype _mD0       = 1.8645;
const fptype piPlusMass = 0.13957018;
const fptype KmMass     = .493677;

int main(int argc, char **argv) {
    GooFit::Application app("Time dependent Dalitz plot, 4 particles", argc, argv);

    TString output = "test_10_15.output";
    app.add_option("-o,--output,output", output, "File to output", true)->check(GooFit::NonexistentPath);

    int trials = 100;
    app.add_option("-t,--trials,output", trials, "Number of trials", true);

    GOOFIT_PARSE(app);

    DecayInfo4t DK3P_DI{Variable("tau", 0.4101, 0.001, 0.300, 0.500),
                        Variable("xmixing", 0.005, 0.001, 0, 0),
                        Variable("ymixing", 0.01, 0.001, 0, 0),
                        Variable("SqWStoRSrate", 1.0 / sqrt(300.0))};

    DK3P_DI.meson_radius = 1.5;
    DK3P_DI.particle_masses.push_back(_mD0);
    DK3P_DI.particle_masses.push_back(piPlusMass);
    DK3P_DI.particle_masses.push_back(piPlusMass);
    DK3P_DI.particle_masses.push_back(KmMass);
    DK3P_DI.particle_masses.push_back(piPlusMass);

    Variable RhoMass{"rho_mass", 0.77526, 0.01, 0.7, 0.8};
    Variable RhoWidth{"rho_width", 0.1478, 0.01, 0.1, 0.2};
    Variable KstarM{"KstarM", 0.89581, 0.01, 0.9, 0.1};
    Variable KstarW{"KstarW", 0.0474, 0.01, 0.1, 0.2};

    // Variable* f600M  = new Variable("f600M", 0.519, 0.01, 0.75, 0.85);
    // Variable* f600W  = new Variable("f600W", 0.454, 0.01, 0.75, 0.85);
    // Variable* a1M  = new Variable("a1M", 1.23, 0.01, 1.2, 1.3);
    // Variable* a1W  = new Variable("a1W", 0.42, 0.01, 0.37, 0.47);
    // Variable* K1M  = new Variable("K1M", 1.272, 0.01, 1.2, 1.3);
    // Variable* K1W  = new Variable("K1W", 0.09, 0.01, 0.08, 0.1);
    // Variable* K1430M  = new Variable("K1430M", 1.414, 0.01, 1.4, 1.5);
    // Variable* K1430W  = new Variable("K1430W", .29, 0.01, 0.25, 0.35);

    // Spin factors: we have two due to the bose symmetrization of the two pi+
    std::vector<SpinFactor *> SFKRS = {new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S, 0, 1, 2, 3),
                                       new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S, 3, 1, 2, 0)};

    std::vector<SpinFactor *> SFKRP;
    SFKRP.push_back(new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_P, 0, 1, 2, 3));
    SFKRP.push_back(new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_P, 3, 1, 2, 0));

    std::vector<SpinFactor *> SFKRD;
    SFKRD.push_back(new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_D, 0, 1, 2, 3));
    SFKRD.push_back(new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_D, 3, 1, 2, 0));

    std::vector<SpinFactor *> SFKF;
    SFKF.push_back(new SpinFactor("SF", SF_4Body::DtoVS_VtoP1P2_StoP3P4, 2, 3, 0, 1));
    SFKF.push_back(new SpinFactor("SF", SF_4Body::DtoVS_VtoP1P2_StoP3P4, 2, 0, 3, 1));

    std::vector<SpinFactor *> SFKK;
    SFKK.push_back(new SpinFactor("SF", SF_4Body::DtoAP1_AtoSP2_StoP3P4, 0, 1, 3, 2));
    SFKK.push_back(new SpinFactor("SF", SF_4Body::DtoAP1_AtoSP2_StoP3P4, 3, 1, 0, 2));

    std::vector<SpinFactor *> SFK1R;
    SFK1R.push_back(new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 3, 2, 0, 1));
    SFK1R.push_back(new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 0, 2, 3, 1));

    std::vector<SpinFactor *> SFA1R;
    SFA1R.push_back(new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 2, 3, 0, 1));
    SFA1R.push_back(new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 2, 0, 3, 1));

    std::vector<SpinFactor *> SFA1RD;
    SFA1RD.push_back(new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2Dwave_VtoP3P4, 2, 3, 0, 1));
    SFA1RD.push_back(new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2Dwave_VtoP3P4, 2, 0, 3, 1));

    // Lineshapes, also for both pi+ configurations
    std::vector<Lineshape *> LSKRS = {new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_12, FF::BL2),
                                      new Lineshapes::RBW("K*(892)bar", KstarM, KstarW, 1, M_34, FF::BL2),
                                      new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_24, FF::BL2),
                                      new Lineshapes::RBW("K*(892)bar", KstarM, KstarW, 1, M_13, FF::BL2)};

    std::vector<Lineshape *> LSKRP = {new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_12, FF::BL2),
                                      new Lineshapes::RBW("K*(892)bar", KstarM, KstarW, 1, M_34, FF::BL2),
                                      new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_24, FF::BL2),
                                      new Lineshapes::RBW("K*(892)bar", KstarM, KstarW, 1, M_13, FF::BL2)};

    std::vector<Lineshape *> LSKRD = {new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_12, FF::BL2),
                                      new Lineshapes::RBW("K*(892)bar", KstarM, KstarW, 1, M_34, FF::BL2),
                                      new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_24, FF::BL2),
                                      new Lineshapes::RBW("K*(892)bar", KstarM, KstarW, 1, M_13, FF::BL2)};

    // the very last parameter means that we have two permutations. so the first half of the Lineshapes
    // and the first half of the spinfactors are amplitude 1, rest is amplitude 2
    // This means that it is important for symmetrized amplitudes that the spinfactors and lineshapes are in the "right"
    // order

    Amplitude Bose_symmetrized_AMP_S{
        "K*(892)rho(770)_S", Variable("amp_real1", 1.0), Variable("amp_imag1", 0.0), LSKRS, SFKRS, 2};
    Amplitude Bose_symmetrized_AMP_P{
        "K*(892)rho(770)_P", Variable("amp_real2", 0.526), Variable("amp_imag2", -0.626), LSKRP, SFKRP, 2};
    Amplitude Bose_symmetrized_AMP_D{
        "K*(892)rho(770)_D", Variable("amp_real3", 26.537), Variable("amp_imag3", 12.284), LSKRD, SFKRD, 2};

    Amplitude Bose_symmetrized_AMP_S_B{
        "B_K*(892)rho(770)_S", Variable("amp_real1", 1.0), Variable("amp_imag1", 0), LSKRS, SFKRS, 2};
    Amplitude Bose_symmetrized_AMP_P_B{
        "B_K*(892)rho(770)_P", Variable("amp_real2", -0.145), Variable("amp_imag2", 0.86), LSKRP, SFKRP, 2};
    Amplitude Bose_symmetrized_AMP_D_B{
        "B_K*(892)rho(770)_D", Variable("amp_real3", 24.343), Variable("amp_imag3", 5.329), LSKRD, SFKRD, 2};

    DK3P_DI.amplitudes_B.push_back(&Bose_symmetrized_AMP_S);
    DK3P_DI.amplitudes_B.push_back(&Bose_symmetrized_AMP_P);
    DK3P_DI.amplitudes_B.push_back(&Bose_symmetrized_AMP_D);

    DK3P_DI.amplitudes.push_back(&Bose_symmetrized_AMP_S_B);
    DK3P_DI.amplitudes.push_back(&Bose_symmetrized_AMP_P_B);
    DK3P_DI.amplitudes.push_back(&Bose_symmetrized_AMP_D_B);

    Observable m12{"m12", 0, 3};
    Observable m34{"m34", 0, 3};
    Observable cos12{"cos12", -1, 1};
    Observable cos34{"m12", -1, 1};
    Observable phi{"phi", -3.5, 3.5};
    EventNumber eventNumber{"eventNumber"};
    Observable dtime{"dtime", 0, 10};
    Observable sigmat{"sigmat", -3, 3};
    Variable constantOne{"constantOne", 1};
    Variable constantZero{"constantZero", 0};

    vector<Observable> observables{m12, m34, cos12, cos34, phi, eventNumber, dtime, sigmat};
    vector<Variable> offsets{constantZero, constantZero};
    vector<Variable> coefficients{constantOne};

    TruthResolution dat;
    PolynomialPdf eff{"constantEff", observables, coefficients, offsets, 0};
    TDDP4 dp{"test", observables, DK3P_DI, &dat, &eff, 0, 1};

    TFile *file = new TFile(output, "RECREATE");
    TTree *tree = new TTree("events", "events");

    double tm12, tm34, tc12, tc34, tphi, tdtime, D0_E, D0_Px, D0_Py, D0_Pz, Kplus_E, Kplus_Px, Kplus_Py, Kplus_Pz,
        Piminus1_E, Piminus1_Px, Piminus1_Py, Piminus1_Pz, Piminus2_E, Piminus2_Px, Piminus2_Py, Piminus2_Pz, Piplus_E,
        Piplus_Px, Piplus_Py, Piplus_Pz;
    int D0_pdg, Kplus_pdg, Piminus1_pdg, Piminus2_pdg, Piplus_pdg;

    tree->Branch("m12", &tm12, "m12/D");
    tree->Branch("m34", &tm34, "m34/D");
    tree->Branch("c12", &tc12, "c12/D");
    tree->Branch("c34", &tc34, "c34/D");
    tree->Branch("phi", &tphi, "phi/D");
    tree->Branch("dtime", &tdtime, "dtime/D");
    tree->Branch("D0_E", &D0_E, "D0_E/D");
    tree->Branch("D0_Px", &D0_Px, "D0_Px/D");
    tree->Branch("D0_Py", &D0_Py, "D0_Py/D");
    tree->Branch("D0_Pz", &D0_Pz, "D0_Pz/D");
    tree->Branch("D0_pdg", &D0_pdg, "D0_pdg/I");
    tree->Branch("Kplus_E", &Kplus_E, "Kplus_E/D");
    tree->Branch("Kplus_Px", &Kplus_Px, "Kplus_Px/D");
    tree->Branch("Kplus_Py", &Kplus_Py, "Kplus_Py/D");
    tree->Branch("Kplus_Pz", &Kplus_Pz, "Kplus_Pz/D");
    tree->Branch("Kplus_pdg", &Kplus_pdg, "Kplus_pdg/I");
    tree->Branch("Piminus1_E", &Piminus1_E, "Piminus1_E/D");
    tree->Branch("Piminus1_Px", &Piminus1_Px, "Piminus1_Px/D");
    tree->Branch("Piminus1_Py", &Piminus1_Py, "Piminus1_Py/D");
    tree->Branch("Piminus1_Pz", &Piminus1_Pz, "Piminus1_Pz/D");
    tree->Branch("Piminus1_pdg", &Piminus1_pdg, "Piminus1_pdg/I");
    tree->Branch("Piminus2_E", &Piminus2_E, "Piminus2_E/D");
    tree->Branch("Piminus2_Px", &Piminus2_Px, "Piminus2_Px/D");
    tree->Branch("Piminus2_Py", &Piminus2_Py, "Piminus2_Py/D");
    tree->Branch("Piminus2_Pz", &Piminus2_Pz, "Piminus2_Pz/D");
    tree->Branch("Piminus2_pdg", &Piminus2_pdg, "Piminus2_pdg/I");
    tree->Branch("Piplus_E", &Piplus_E, "Piplus_E/D");
    tree->Branch("Piplus_Px", &Piplus_Px, "Piplus_Px/D");
    tree->Branch("Piplus_Py", &Piplus_Py, "Piplus_Py/D");
    tree->Branch("Piplus_Pz", &Piplus_Pz, "Piplus_Pz/D");
    tree->Branch("Piplus_pdg", &Piplus_pdg, "Piplus_pdg/I");

    int total_accepted = 0;

    for(int k = 0; k < trials; ++k) {
        int numEvents = 800000;
        dp.setGenerationOffset(k * numEvents);

        mcbooster::ParticlesSet_h particles; // typedef for std::vector<Particles_h *>
        mcbooster::VariableSet_h variables;
        mcbooster::RealVector_h weights;
        mcbooster::BoolVector_h flags;

        std::tie(particles, variables, weights, flags) = dp.GenerateSig(numEvents);

        int accepted = thrust::count_if(flags.begin(), flags.end(), thrust::identity<bool>());
        total_accepted += accepted;

        GOOFIT_INFO(
            "Run #{}: Using accept-reject method would leave you with {} out of {} events", k, accepted, numEvents);

        for(int i = 0; i < weights.size(); ++i) {
            if(flags[i] == 1) {
                // printf("%.5g %.5g %.5g %.5g %.5g %.5g %.5g\n", (*(variables[0]))[i], (*(variables[1]))[i],
                // (*(variables[2]))[i], (*(variables[3]))[i], (*(variables[4]))[i], weights[i], flags[i]);
                tm12         = (*(variables[0]))[i];
                tm34         = (*(variables[1]))[i];
                tc12         = (*(variables[2]))[i];
                tc34         = (*(variables[3]))[i];
                tphi         = (*(variables[4]))[i];
                tdtime       = (*(variables[5]))[i];
                D0_E         = 1864;
                D0_Px        = 0.0;
                D0_Py        = 0.0;
                D0_Pz        = 0.0;
                D0_pdg       = 421;
                Kplus_E      = 1000 * (*(particles[2]))[i].get(0);
                Kplus_Px     = 1000 * (*(particles[2]))[i].get(1);
                Kplus_Py     = 1000 * (*(particles[2]))[i].get(2);
                Kplus_Pz     = 1000 * (*(particles[2]))[i].get(3);
                Kplus_pdg    = -321;
                Piminus1_E   = 1000 * (*(particles[3]))[i].get(0);
                Piminus1_Px  = 1000 * (*(particles[3]))[i].get(1);
                Piminus1_Py  = 1000 * (*(particles[3]))[i].get(2);
                Piminus1_Pz  = 1000 * (*(particles[3]))[i].get(3);
                Piminus1_pdg = 211;
                Piminus2_E   = 1000 * (*(particles[0]))[i].get(0);
                Piminus2_Px  = 1000 * (*(particles[0]))[i].get(1);
                Piminus2_Py  = 1000 * (*(particles[0]))[i].get(2);
                Piminus2_Pz  = 1000 * (*(particles[0]))[i].get(3);
                Piminus2_pdg = 211;
                Piplus_E     = 1000 * (*(particles[1]))[i].get(0);
                Piplus_Px    = 1000 * (*(particles[1]))[i].get(1);
                Piplus_Py    = 1000 * (*(particles[1]))[i].get(2);
                Piplus_Pz    = 1000 * (*(particles[1]))[i].get(3);
                Piplus_pdg   = -211;

                tree->Fill();
            }
        }

        delete variables[0];
        delete variables[1];
        delete variables[2];
        delete variables[3];
        delete variables[4];
        delete variables[5];

        delete particles[0];
        delete particles[1];
        delete particles[2];
        delete particles[3];
    }

    tree->Write();
    file->Close();

    if(total_accepted > 0)
        return 0;
    else {
        GOOFIT_ERROR("Total accepted was 0! Something is wrong.");
        return 1;
    }
}
