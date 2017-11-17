#include <fstream>

// GooFit stuff
#include <goofit/Application.h>
#include <goofit/FitManager.h>
#include <goofit/PDFs/basic/PolynomialPdf.h>
#include <goofit/PDFs/combine/AddPdf.h>
#include <goofit/PDFs/physics/DP4Pdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>

using namespace std;
using namespace GooFit;

const fptype _mD0       = 1.8645;
const fptype piPlusMass = 0.13957018;
const fptype KmMass     = .493677;
// Constants used in more than one PDF component.

int main(int argc, char **argv) {
    GooFit::Application app("Dalitz 4 daughter example", argc, argv);

    // Set all host normalisation just to make sure no errors
    for(int i = 0; i < maxParams; i++)
        host_normalisation[i] = -7;

    GOOFIT_PARSE(app, argc, argv);

    Observable m12{"m12", 0, 3};
    Observable m34{"m34", 0, 3};
    Observable cos12{"cos12", -1, 1};
    Observable cos34{"m12", -1, 1};
    Observable phi{"phi", -3.5, 3.5};
    EventNumber eventNumber{"eventNumber", 0, INT_MAX};

    UnbinnedDataSet currData{m12, m34, cos12, cos34, phi, eventNumber};

    unsigned int MCevents = 0;

    std::string input_str = app.get_filename("ToyMC.txt", "examples/DP4");

    { // FStream block
        fstream input(input_str, std::ios_base::in);

        while(input >> m12 >> m34 >> cos12 >> cos34 >> phi) {
            eventNumber.setValue(MCevents++);
            currData.addEvent();
        }
    }

    GOOFIT_INFO("Read in {} events", MCevents);

    DecayInfo4 DK3P_DI;
    DK3P_DI.meson_radius    = 1.5;
    DK3P_DI.particle_masses = {_mD0, piPlusMass, piPlusMass, KmMass, piPlusMass};

    Variable RhoMass{"rho_mass", 0.77526};
    Variable RhoWidth{"rho_width", 0.1478};
    Variable KstarM{"KstarM", 0.89581};
    Variable KstarW{"KstarW", 0.0474};
    Variable f600M{"f600M", 0.519};
    Variable f600W{"f600W", 0.454};
    Variable a1M{"a1M", 1.23};
    Variable a1W{"a1W", 0.42};
    Variable K1M{"K1M", 1.272};
    Variable K1W{"K1W", 0.09};
    Variable K1430M{"K1430M", 1.414};
    Variable K1430W{"K1430W", 0.29};

    // Spin factors: we have two due to the bose symmetrization of the two pi+
    std::vector<SpinFactor *> SFKRS = {new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S, 0, 1, 2, 3),
                                       new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S, 3, 1, 2, 0)};

    std::vector<SpinFactor *> SFKRP = {new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_P, 0, 1, 2, 3),
                                       new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_P, 3, 1, 2, 0)};

    std::vector<SpinFactor *> SFKRD = {new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_D, 0, 1, 2, 3),
                                       new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_D, 3, 1, 2, 0)};

    std::vector<SpinFactor *> SFKF = {new SpinFactor("SF", SF_4Body::DtoVS_VtoP1P2_StoP3P4, 2, 3, 0, 1),
                                      new SpinFactor("SF", SF_4Body::DtoVS_VtoP1P2_StoP3P4, 2, 0, 3, 1)};

    std::vector<SpinFactor *> SFKK = {new SpinFactor("SF", SF_4Body::DtoAP1_AtoSP2_StoP3P4, 0, 1, 3, 2),
                                      new SpinFactor("SF", SF_4Body::DtoAP1_AtoSP2_StoP3P4, 3, 1, 0, 2)};

    std::vector<SpinFactor *> SFK1R = {new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 3, 2, 0, 1),
                                       new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 0, 2, 3, 1)};

    std::vector<SpinFactor *> SFA1R = {new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 2, 3, 0, 1),
                                       new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 2, 0, 3, 1)};

    std::vector<SpinFactor *> SFA1RD = {new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2Dwave_VtoP3P4, 2, 3, 0, 1),
                                        new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2Dwave_VtoP3P4, 2, 0, 3, 1)};

    // Lineshapes, also for both pi+ configurations
    std::vector<Lineshape *> LSKRS = {new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_12),
                                      new Lineshapes::RBW("K*(892)bar", KstarM, KstarW, 1, M_34),
                                      new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_24),
                                      new Lineshapes::RBW("K*(892)bar", KstarM, KstarW, 1, M_13)};

    std::vector<Lineshape *> LSKRP = {new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_12),
                                      new Lineshapes::RBW("K*(892)bar", KstarM, KstarW, 1, M_34),
                                      new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_24),
                                      new Lineshapes::RBW("K*(892)bar", KstarM, KstarW, 1, M_13)};

    std::vector<Lineshape *> LSKRD = {new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_12),
                                      new Lineshapes::RBW("K*(892)bar", KstarM, KstarW, 1, M_34),
                                      new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_24),
                                      new Lineshapes::RBW("K*(892)bar", KstarM, KstarW, 1, M_13)};

    std::vector<Lineshape *> LSKF = {new Lineshapes::RBW("K*(892)bar", KstarM, KstarW, 1, M_34),
                                     new Lineshapes::Bugg("f600", f600M, f600W, 0, M_12),
                                     new Lineshapes::RBW("K*(892)bar", KstarM, KstarW, 1, M_13),
                                     new Lineshapes::Bugg("f600", f600M, f600W, 0, M_24)};

    std::vector<Lineshape *> LSKK = {new Lineshapes::SBW("K(1)(1270)bar", K1M, K1W, 1, M_34_2),
                                     new Lineshapes::LASS("K(0)*(1430)bar", K1430M, K1430W, 0, M_34),
                                     new Lineshapes::SBW("K(1)(1270)bar2", K1M, K1W, 1, M_13_2),
                                     new Lineshapes::LASS("K(0)*(1430)bar2", K1430M, K1430W, 0, M_13)};

    std::vector<Lineshape *> LSK1R = {new Lineshapes::SBW("K(1)(1270)bar", K1M, K1W, 0, M_12_3),
                                      new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_12),
                                      new Lineshapes::SBW("K(1)(1270)bar", K1M, K1W, 0, M_24_3),
                                      new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_24)};

    std::vector<Lineshape *> LSA1R = {new Lineshapes::SBW("a(1)(1260)+", a1M, a1W, 0, M_12_4),
                                      new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_12),
                                      new Lineshapes::SBW("a(1)(1260)+", a1M, a1W, 0, M_24_1),
                                      new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_24)};

    std::vector<Lineshape *> LSA1RD = {new Lineshapes::SBW("a(1)(1260)+", a1M, a1W, 2, M_12_4),
                                       new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_12),
                                       new Lineshapes::SBW("a(1)(1260)+", a1M, a1W, 2, M_24_1),
                                       new Lineshapes::RBW("rho(770)", RhoMass, RhoWidth, 1, M_24)};

    // the very last parameter means that we have two permutations. so the first half of the Lineshapes
    // and the first half of the spinfactors are amplitude 1, rest is amplitude
    // This means that it is important for symmetrized amplitueds that the spinfactors and lineshapes are in the "right"
    // order

    // Amplitudes with the correct fixed values.
    // Amplitude* Bose_symmetrized_AMP_S = new Amplitude( "K*(892)rho(770)_S", new Variable("amp_real1", -0.115177), new
    // Variable("amp_imag1", 0.153976), LSKRS, SFKRS, 2);
    // Amplitude* Bose_symmetrized_AMP_P = new Amplitude( "K*(892)rho(770)_P", new Variable("amp_real2", -0.0298697),
    // new Variable("amp_imag2", -0.0722874), LSKRP, SFKRP, 2);
    // Amplitude* Bose_symmetrized_AMP_D = new Amplitude( "K*(892)rho(770)_D", new Variable("amp_real3", -0.452212), new
    // Variable("amp_imag3", 0.426521), LSKRD, SFKRD, 2);

    // Amplitudes with floating slightly different values to be fitted.

    Amplitude Bose_symmetrized_AMP_S{"K*(892)rho(770)_S",
                                     Variable("amp_real1", -0.1, 0.001, 0, 0),
                                     Variable("amp_imag1", 0.1, 0.001, 0, 0),
                                     LSKRS,
                                     SFKRS,
                                     2};

    Amplitude Bose_symmetrized_AMP_P{"K*(892)rho(770)_P",
                                     Variable("amp_real2", -0.02, 0.001, 0, 0),
                                     Variable("amp_imag2", -0.07, 0.001, 0, 0),
                                     LSKRP,
                                     SFKRP,
                                     2};
    Amplitude Bose_symmetrized_AMP_D{"K*(892)rho(770)_D",
                                     Variable("amp_real3", -0.4, 0.001, 0, 0),
                                     Variable("amp_imag3", 0.4, 0.001, 0, 0),
                                     LSKRD,
                                     SFKRD,
                                     2};

    Amplitude Bose_symmetrized_KF{
        "KF", Variable("amp_real4", 0.0120787), Variable("amp_imag4", -0.0332525), LSKF, SFKF, 2};
    Amplitude Bose_symmetrized_KK{
        "LSKK", Variable("amp_real5", 0.0109033), Variable("amp_imag5", -0.00186219), LSKK, SFKK, 2};
    Amplitude Bose_symmetrized_K1R{
        "LSK1R", Variable("amp_real6", -0.10728), Variable("amp_imag6", -0.130213), LSK1R, SFK1R, 2};
    Amplitude Bose_symmetrized_A1R{"LSA1R", Variable("amp_real7", 1.0), Variable("amp_imag7", 0.0), LSA1R, SFA1R, 2};
    Amplitude Bose_symmetrized_A1RD{
        "LSA1RD", Variable("amp_real8", -0.94921), Variable("amp_imag8", -1.73407), LSA1RD, SFA1RD, 2};

    DK3P_DI.amplitudes = {&Bose_symmetrized_KF,
                          &Bose_symmetrized_AMP_S,
                          &Bose_symmetrized_AMP_P,
                          &Bose_symmetrized_AMP_D,
                          &Bose_symmetrized_KK,
                          &Bose_symmetrized_K1R,
                          &Bose_symmetrized_A1R,
                          &Bose_symmetrized_A1RD};

    for(auto res : LSKF)
        res->setParameterConstantness(true);

    for(auto res : LSKRS)
        res->setParameterConstantness(true);

    for(auto res : LSKRP)
        res->setParameterConstantness(true);

    for(auto res : LSKRD)
        res->setParameterConstantness(true);

    for(auto res : LSKK)
        res->setParameterConstantness(true);

    for(auto res : LSK1R)
        res->setParameterConstantness(true);

    for(auto res : LSA1R)
        res->setParameterConstantness(true);

    for(auto res : LSA1RD)
        res->setParameterConstantness(true);

    Variable constantOne{"constantOne", 1};
    Variable constantZero{"constantZero", 0};

    vector<Observable> observables = {m12, m34, cos12, cos34, phi, eventNumber};
    vector<Variable> coefficients  = {constantOne};
    vector<Variable> offsets       = {constantZero, constantZero};

    PolynomialPdf eff{"constantEff", observables, coefficients, offsets, 0};
    DPPdf dp{"test", observables, DK3P_DI, &eff, 1000000};

    Variable constant{"constant", 0.1};
    Variable constant2{"constant2", 1.0};
    PolynomialPdf backgr{"backgr", m12, {constant}};
    AddPdf signal{"signal", constant2, &dp, &backgr};

    signal.setData(&currData);

    dp.setDataSize(currData.getNumEvents(), 6);

    FitManager datapdf{&signal};
    datapdf.fit();

    return datapdf;
}
