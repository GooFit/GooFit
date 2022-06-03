// ROOT
#include <TFile.h>
#include <TTree.h>

// GooFit stuff
#include <goofit/Application.h>
#include <goofit/PDFs/basic/PolynomialPdf.h>
#include <goofit/PDFs/combine/AddPdf.h>
#include <goofit/PDFs/physics/DP4Pdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>

#include <thrust/count.h>

using namespace GooFit;

// Constants used in more than one PDF component.
const fptype _mD0       = 1.8645;
const fptype piPlusMass = 0.13957018;
const fptype KmMass     = .493677;

int main(int argc, char **argv) {
    GooFit::Application app("Signal Generator Example", argc, argv);

    GOOFIT_PARSE(app);

    DecayInfo4 DK3P_DI;
    DK3P_DI.meson_radius = 1.5;
    DK3P_DI.particle_masses.push_back(_mD0);
    DK3P_DI.particle_masses.push_back(piPlusMass);
    DK3P_DI.particle_masses.push_back(piPlusMass);
    DK3P_DI.particle_masses.push_back(KmMass);
    DK3P_DI.particle_masses.push_back(piPlusMass);

    Variable RhoMass("rho_mass", 0.77526, 0.01, 0.7, 0.8);
    Variable RhoWidth("rho_width", 0.1478, 0.01, 0.1, 0.2);
    Variable KstarM("KstarM", 0.89581, 0.01, 0.9, 0.1);
    Variable KstarW("KstarW", 0.0474, 0.01, 0.1, 0.2);
    Variable f600M("f600M", 0.519, 0.01, 0.75, 0.85);
    Variable f600W("f600W", 0.454, 0.01, 0.75, 0.85);
    Variable a1M("a1M", 1.23, 0.01, 1.2, 1.3);
    Variable a1W("a1W", 0.42, 0.01, 0.37, 0.47);
    Variable K1M("K1M", 1.272, 0.01, 1.2, 1.3);
    Variable K1W("K1W", 0.09, 0.01, 0.08, 0.1);
    Variable K1430M("K1430M", 1.414, 0.01, 1.4, 1.5);
    Variable K1430W("K1430W", .29, 0.01, 0.25, 0.35);

    // Spin factors: we have two due to the bose symmetrization of the two pi+
    std::vector<SpinFactor *> SFKRS;
    SFKRS.push_back(new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S, _mD0, 0, 1, 2, 3));
    SFKRS.push_back(new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S, _mD0, 3, 1, 2, 0));

    std::vector<SpinFactor *> SFKRP;
    SFKRP.push_back(new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_P, _mD0, 0, 1, 2, 3));
    SFKRP.push_back(new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_P, _mD0, 3, 1, 2, 0));

    std::vector<SpinFactor *> SFKRD;
    SFKRD.push_back(new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_D, _mD0, 0, 1, 2, 3));
    SFKRD.push_back(new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_D, _mD0, 3, 1, 2, 0));

    std::vector<SpinFactor *> SFKF;
    SFKF.push_back(new SpinFactor("SF", SF_4Body::DtoVS_VtoP1P2_StoP3P4, _mD0, 2, 3, 0, 1));
    SFKF.push_back(new SpinFactor("SF", SF_4Body::DtoVS_VtoP1P2_StoP3P4, _mD0, 2, 0, 3, 1));

    std::vector<SpinFactor *> SFKK;
    SFKK.push_back(new SpinFactor("SF", SF_4Body::DtoAP1_AtoSP2_StoP3P4, _mD0, 0, 1, 3, 2));
    SFKK.push_back(new SpinFactor("SF", SF_4Body::DtoAP1_AtoSP2_StoP3P4, _mD0, 3, 1, 0, 2));

    std::vector<SpinFactor *> SFK1R;
    SFK1R.push_back(new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, _mD0, 3, 2, 0, 1));
    SFK1R.push_back(new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, _mD0, 0, 2, 3, 1));

    std::vector<SpinFactor *> SFA1R;
    SFA1R.push_back(new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, _mD0, 2, 3, 0, 1));
    SFA1R.push_back(new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, _mD0, 2, 0, 3, 1));

    std::vector<SpinFactor *> SFA1RD;
    SFA1RD.push_back(new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2Dwave_VtoP3P4, _mD0, 2, 3, 0, 1));
    SFA1RD.push_back(new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2Dwave_VtoP3P4, _mD0, 2, 0, 3, 1));

    // Lineshapes, also for both pi+ configurations
    std::vector<Lineshape *> LSKRS;
    LSKRS.push_back(new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_12));
    LSKRS.push_back(new Lineshapes::RBW("K*(892)bar", KstarM, KstarW, 1, M_34));
    LSKRS.push_back(new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_24));
    LSKRS.push_back(new Lineshapes::RBW("K*(892)bar", KstarM, KstarW, 1, M_13));

    std::vector<Lineshape *> LSKRP;
    LSKRP.push_back(new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_12));
    LSKRP.push_back(new Lineshapes::RBW("K*(892)bar", KstarM, KstarW, 1, M_34));
    LSKRP.push_back(new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_24));
    LSKRP.push_back(new Lineshapes::RBW("K*(892)bar", KstarM, KstarW, 1, M_13));

    std::vector<Lineshape *> LSKRD;
    LSKRD.push_back(new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_12));
    LSKRD.push_back(new Lineshapes::RBW("K*(892)bar", KstarM, KstarW, 1, M_34));
    LSKRD.push_back(new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_24));
    LSKRD.push_back(new Lineshapes::RBW("K*(892)bar", KstarM, KstarW, 1, M_13));

    std::vector<Lineshape *> LSKF;
    LSKF.push_back(new Lineshapes::RBW("K*(892)bar", KstarM, KstarW, 1, M_34));
    LSKF.push_back(new Lineshapes::Bugg("f600", f600M, f600W, 0, M_12));
    LSKF.push_back(new Lineshapes::RBW("K*(892)bar", KstarM, KstarW, 1, M_13));
    LSKF.push_back(new Lineshapes::Bugg("f600", f600M, f600W, 0, M_24));

    std::vector<Lineshape *> LSKK;
    LSKK.push_back(new Lineshapes::SBW("K(1)(1270)bar", K1M, K1W, 1, M_34_2));
    LSKK.push_back(new Lineshapes::LASS("K(0)*(1430)bar", K1430M, K1430W, 0, M_34));
    LSKK.push_back(new Lineshapes::SBW("K(1)(1270)bar2", K1M, K1W, 1, M_13_2));
    LSKK.push_back(new Lineshapes::LASS("K(0)*(1430)bar2", K1430M, K1430W, 0, M_13));

    std::vector<Lineshape *> LSK1R;
    LSK1R.push_back(new Lineshapes::SBW("K(1)(1270)bar", K1M, K1W, 0, M_12_3));
    LSK1R.push_back(new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_12));
    LSK1R.push_back(new Lineshapes::SBW("K(1)(1270)bar", K1M, K1W, 0, M_24_3));
    LSK1R.push_back(new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_24));

    std::vector<Lineshape *> LSA1R;
    LSA1R.push_back(new Lineshapes::SBW("a(1)(1260)+", a1M, a1W, 0, M_12_4));
    LSA1R.push_back(new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_12));
    LSA1R.push_back(new Lineshapes::SBW("a(1)(1260)+", a1M, a1W, 0, M_24_1));
    LSA1R.push_back(new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_24));

    std::vector<Lineshape *> LSA1RD;
    LSA1RD.push_back(new Lineshapes::SBW("a(1)(1260)+", a1M, a1W, 2, M_12_4));
    LSA1RD.push_back(new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_12));
    LSA1RD.push_back(new Lineshapes::SBW("a(1)(1260)+", a1M, a1W, 2, M_24_1));
    LSA1RD.push_back(new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_24));

    // the very last parameter means that we have two permutations. so the first half of the Lineshapes
    // and the first half of the spinfactors are amplitude 1, rest is amplitude 2
    // This means that it is important for symmetrized amplitueds that the spinfactors and lineshapes are in the "right"
    // order
    Amplitude *Bose_symmetrized_AMP_S = new Amplitude(
        "K*(892)rho(770)_S", Variable("amp_real1", -0.115177), Variable("amp_imag1", 0.153976), LSKRS, SFKRS, 2);
    Amplitude *Bose_symmetrized_AMP_P = new Amplitude(
        "K*(892)rho(770)_P", Variable("amp_real2", -0.0298697), Variable("amp_imag2", -0.0722874), LSKRP, SFKRP, 2);
    Amplitude *Bose_symmetrized_AMP_D = new Amplitude(
        "K*(892)rho(770)_D", Variable("amp_real3", -0.452212), Variable("amp_imag3", 0.426521), LSKRD, SFKRD, 2);

    Amplitude *Bose_symmetrized_KF
        = new Amplitude("KF", Variable("amp_real3", 0.0120787), Variable("amp_imag3", -0.0332525), LSKF, SFKF, 2);
    Amplitude *Bose_symmetrized_KK
        = new Amplitude("LSKK", Variable("amp_real3", 0.0109033), Variable("amp_imag3", -0.00186219), LSKK, SFKK, 2);
    Amplitude *Bose_symmetrized_K1R
        = new Amplitude("LSK1R", Variable("amp_real3", -0.10728), Variable("amp_imag3", -0.130213), LSK1R, SFK1R, 2);
    Amplitude *Bose_symmetrized_A1R
        = new Amplitude("LSA1R", Variable("amp_real3", 1.0), Variable("amp_imag3", 0.0), LSA1R, SFA1R, 2);
    Amplitude *Bose_symmetrized_A1RD
        = new Amplitude("LSA1RD", Variable("amp_real3", -0.94921), Variable("amp_imag3", -1.73407), LSA1RD, SFA1RD, 2);

    DK3P_DI.amplitudes.push_back(Bose_symmetrized_KF);
    DK3P_DI.amplitudes.push_back(Bose_symmetrized_AMP_S);
    DK3P_DI.amplitudes.push_back(Bose_symmetrized_AMP_P);
    DK3P_DI.amplitudes.push_back(Bose_symmetrized_AMP_D);
    DK3P_DI.amplitudes.push_back(Bose_symmetrized_KK);
    DK3P_DI.amplitudes.push_back(Bose_symmetrized_K1R);
    DK3P_DI.amplitudes.push_back(Bose_symmetrized_A1R);
    DK3P_DI.amplitudes.push_back(Bose_symmetrized_A1RD);

    Observable m12("m12", 0, 3);
    Observable m34("m34", 0, 3);
    Observable cos12("cos12", -1, 1);
    Observable cos34("m12", -1, 1);
    Observable phi("phi", -3.5, 3.5);
    EventNumber eventNumber("eventNumber");
    Variable constantOne("constantOne", 1);
    Variable constantZero("constantZero", 0);

    std::vector<Observable> observables;
    std::vector<Variable> coefficients;
    std::vector<Variable> offsets;

    observables.push_back(m12);
    observables.push_back(m34);
    observables.push_back(cos12);
    observables.push_back(cos34);
    observables.push_back(phi);
    observables.push_back(eventNumber);
    offsets.push_back(constantZero);
    offsets.push_back(constantZero);
    coefficients.push_back(constantOne);

    PolynomialPdf eff("constantEff", observables, coefficients, offsets, 0);
    auto dp = new DPPdf("test", observables, DK3P_DI, &eff, 5);

    TFile file("SigGen.root", "RECREATE");
    TTree tree("events", "events");

    double tm12, tm34, tc12, tc34, tphi;
    tree.Branch("m12", &tm12, "m12/D");
    tree.Branch("m34", &tm34, "m34/D");
    tree.Branch("c12", &tc12, "c12/D");
    tree.Branch("c34", &tc34, "c34/D");
    tree.Branch("phi", &tphi, "phi/D");

    for(int k = 0; k < 4; ++k) {
        int numEvents = 1e6;
        dp->setGenerationOffset(k * numEvents);
        auto tuple = dp->GenerateSig(numEvents);

        auto particles = std::get<0>(tuple);
        auto variables = std::get<1>(tuple);
        auto weights   = std::get<2>(tuple);
        auto flags     = std::get<3>(tuple);
        long accepted  = thrust::count_if(flags.begin(), flags.end(), thrust::identity<bool>());
        fmt::print(
            stderr, "Using accept-reject method would leave you with {} out of {} events\n", accepted, numEvents);

        for(int i = 0; i < weights.size(); ++i) {
            if(flags[i] == 1) {
                tm12 = (*(variables[0]))[i];
                tm34 = (*(variables[1]))[i];
                tc12 = (*(variables[2]))[i];
                tc34 = (*(variables[3]))[i];
                tphi = (*(variables[4]))[i];
                tree.Fill();
            }
        }

        delete variables[0];
        delete variables[1];
        delete variables[2];
        delete variables[3];
        delete variables[4];

        delete particles[0];
        delete particles[1];
        delete particles[2];
        delete particles[3];

        if(accepted == 0) {
            GOOFIT_ERROR("ERROR: 0 Events accepted! Something is wrong with PDF evaluation (probably)");
            return 2;
        }
    }

    tree.Write();
    return 0;
}
