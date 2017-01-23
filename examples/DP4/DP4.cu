#include <fstream>

// GooFit stuff
#include "goofit/FitManager.h"
#include "goofit/Variable.h"
#include "goofit/PDFs/PolynomialPdf.h"
#include "goofit/PDFs/AddPdf.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/PDFs/DP4Pdf.h"
using namespace std;

const fptype _mD0 = 1.8645;
const fptype piPlusMass = 0.13957018;
const fptype KmMass = .493677;
// Constants used in more than one PDF component.

int main(int argc, char** argv) {

    Variable* m12 = new Variable("m12", 0, 3);
    Variable* m34 = new Variable("m34", 0, 3);
    Variable* cos12 = new Variable("cos12", -1, 1);
    Variable* cos34 = new Variable("m12", -1, 1);
    Variable* phi = new Variable("phi", -3.5, 3.5);
    Variable* eventNumber = new Variable("eventNumber", 0, INT_MAX);

    std::vector<Variable*> vars;
    vars.push_back(m12);
    vars.push_back(m34);
    vars.push_back(cos12);
    vars.push_back(cos34);
    vars.push_back(phi);
    vars.push_back(eventNumber);
    UnbinnedDataSet currData(vars);

    unsigned int MCevents = 0;

    fstream input("ToyMC.txt", std::ios_base::in);

    while(input >> m12->value >> m34->value >> cos12->value >> cos34->value >> phi->value) {
        eventNumber->value = MCevents++;
        currData.addEvent();
    }

    printf("read in %i events\n", MCevents);

    DecayInfo_DP* DK3P_DI = new DecayInfo_DP();
    DK3P_DI->meson_radius =1.5;
    DK3P_DI->particle_masses.push_back(_mD0);
    DK3P_DI->particle_masses.push_back(piPlusMass);
    DK3P_DI->particle_masses.push_back(piPlusMass);
    DK3P_DI->particle_masses.push_back(KmMass);
    DK3P_DI->particle_masses.push_back(piPlusMass);

    Variable* RhoMass  = new Variable("rho_mass", 0.77526);
    Variable* RhoWidth = new Variable("rho_width", 0.1478);
    Variable* KstarM   = new Variable("KstarM", 0.89581);
    Variable* KstarW   = new Variable("KstarW", 0.0474);
    Variable* f600M    = new Variable("f600M", 0.519);
    Variable* f600W    = new Variable("f600W", 0.454);
    Variable* a1M      = new Variable("a1M", 1.23);
    Variable* a1W      = new Variable("a1W", 0.42);
    Variable* K1M      = new Variable("K1M", 1.272);
    Variable* K1W      = new Variable("K1W", 0.09);
    Variable* K1430M   = new Variable("K1430M", 1.414);
    Variable* K1430W   = new Variable("K1430W", 0.29);

    //Spin factors: we have two due to the bose symmetrization of the two pi+
    std::vector<SpinFactor*> SFKRS;
    SFKRS.push_back(new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S, 0, 1, 2, 3));
    SFKRS.push_back(new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S, 3, 1, 2, 0));

    std::vector<SpinFactor*> SFKRP;
    SFKRP.push_back(new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_P, 0, 1, 2, 3));
    SFKRP.push_back(new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_P, 3, 1, 2, 0));

    std::vector<SpinFactor*> SFKRD;
    SFKRD.push_back(new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_D, 0, 1, 2, 3));
    SFKRD.push_back(new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_D, 3, 1, 2, 0));

    std::vector<SpinFactor*> SFKF;
    SFKF.push_back(new SpinFactor("SF", SF_4Body::DtoVS_VtoP1P2_StoP3P4, 2, 3, 0, 1));
    SFKF.push_back(new SpinFactor("SF", SF_4Body::DtoVS_VtoP1P2_StoP3P4, 2, 0, 3, 1));

    std::vector<SpinFactor*> SFKK;
    SFKK.push_back(new SpinFactor("SF", SF_4Body::DtoAP1_AtoSP2_StoP3P4, 0, 1, 3, 2));
    SFKK.push_back(new SpinFactor("SF", SF_4Body::DtoAP1_AtoSP2_StoP3P4, 3, 1, 0, 2));

    std::vector<SpinFactor*> SFK1R;
    SFK1R.push_back(new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 3, 2, 0, 1));
    SFK1R.push_back(new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 0, 2, 3, 1));

    std::vector<SpinFactor*> SFA1R;
    SFA1R.push_back(new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 2, 3, 0, 1));
    SFA1R.push_back(new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 2, 0, 3, 1));

    std::vector<SpinFactor*> SFA1RD;
    SFA1RD.push_back(new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2Dwave_VtoP3P4, 2, 3, 0, 1));
    SFA1RD.push_back(new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2Dwave_VtoP3P4, 2, 0, 3, 1));

    //Lineshapes, also for both pi+ configurations
    std::vector<Lineshape*> LSKRS;
    LSKRS.push_back(new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, LS::BW));
    LSKRS.push_back(new Lineshape("K*(892)bar", KstarM, KstarW, 1, M_34, LS::BW));
    LSKRS.push_back(new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, LS::BW));
    LSKRS.push_back(new Lineshape("K*(892)bar", KstarM, KstarW, 1, M_13, LS::BW));

    std::vector<Lineshape*> LSKRP;
    LSKRP.push_back(new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, LS::BW));
    LSKRP.push_back(new Lineshape("K*(892)bar", KstarM, KstarW, 1, M_34, LS::BW));
    LSKRP.push_back(new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, LS::BW));
    LSKRP.push_back(new Lineshape("K*(892)bar", KstarM, KstarW, 1, M_13, LS::BW));

    std::vector<Lineshape*> LSKRD;
    LSKRD.push_back(new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, LS::BW));
    LSKRD.push_back(new Lineshape("K*(892)bar", KstarM, KstarW, 1, M_34, LS::BW));
    LSKRD.push_back(new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, LS::BW));
    LSKRD.push_back(new Lineshape("K*(892)bar", KstarM, KstarW, 1, M_13, LS::BW));

    std::vector<Lineshape*> LSKF;
    LSKF.push_back(new Lineshape("K*(892)bar", KstarM, KstarW, 1, M_34, LS::BW));
    LSKF.push_back(new Lineshape("f600", f600M, f600W, 0, M_12, LS::Bugg));
    LSKF.push_back(new Lineshape("K*(892)bar", KstarM, KstarW, 1, M_13, LS::BW));
    LSKF.push_back(new Lineshape("f600", f600M, f600W, 0, M_24, LS::Bugg));

    std::vector<Lineshape*> LSKK;
    LSKK.push_back(new Lineshape("K(1)(1270)bar", K1M, K1W, 1, M_34_2, LS::SBW));
    LSKK.push_back(new Lineshape("K(0)*(1430)bar", K1430M, K1430W, 0, M_34, LS::Lass));
    LSKK.push_back(new Lineshape("K(1)(1270)bar2", K1M, K1W, 1, M_13_2, LS::SBW));
    LSKK.push_back(new Lineshape("K(0)*(1430)bar2", K1430M, K1430W, 0, M_13, LS::Lass));

    std::vector<Lineshape*> LSK1R;
    LSK1R.push_back(new Lineshape("K(1)(1270)bar", K1M, K1W, 0, M_12_3, LS::SBW));
    LSK1R.push_back(new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, LS::BW));
    LSK1R.push_back(new Lineshape("K(1)(1270)bar", K1M, K1W, 0, M_24_3, LS::SBW));
    LSK1R.push_back(new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, LS::BW));

    std::vector<Lineshape*> LSA1R;
    LSA1R.push_back(new Lineshape("a(1)(1260)+", a1M, a1W, 0, M_12_4, LS::SBW));
    LSA1R.push_back(new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, LS::BW));
    LSA1R.push_back(new Lineshape("a(1)(1260)+", a1M, a1W, 0, M_24_1, LS::SBW));
    LSA1R.push_back(new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, LS::BW));

    std::vector<Lineshape*> LSA1RD;
    LSA1RD.push_back(new Lineshape("a(1)(1260)+", a1M, a1W, 2, M_12_4, LS::SBW));
    LSA1RD.push_back(new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, LS::BW));
    LSA1RD.push_back(new Lineshape("a(1)(1260)+", a1M, a1W, 2, M_24_1, LS::SBW));
    LSA1RD.push_back(new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, LS::BW));


    // the very last parameter means that we have two permutations. so the first half of the Lineshapes
    // and the first half of the spinfactors are amplitude 1, rest is amplitude
    // This means that it is important for symmetrized amplitueds that the spinfactors and lineshapes are in the "right" order

    // Amplitudes with the correct fixed values.
    // Amplitude* Bose_symmetrized_AMP_S = new Amplitude( "K*(892)rho(770)_S", new Variable("amp_real1", -0.115177), new Variable("amp_imag1", 0.153976), LSKRS, SFKRS, 2);
    // Amplitude* Bose_symmetrized_AMP_P = new Amplitude( "K*(892)rho(770)_P", new Variable("amp_real2", -0.0298697), new Variable("amp_imag2", -0.0722874), LSKRP, SFKRP, 2);
    // Amplitude* Bose_symmetrized_AMP_D = new Amplitude( "K*(892)rho(770)_D", new Variable("amp_real3", -0.452212), new Variable("amp_imag3", 0.426521), LSKRD, SFKRD, 2);

    //Amplitudes with floating slightliy different values to be fitted.
    Amplitude* Bose_symmetrized_AMP_S = new Amplitude("K*(892)rho(770)_S", new Variable("amp_real1", -0.1, 0.001, 0, 0),
            new Variable("amp_imag1", 0.1, 0.001, 0, 0), LSKRS, SFKRS, 2);
    Amplitude* Bose_symmetrized_AMP_P = new Amplitude("K*(892)rho(770)_P", new Variable("amp_real2", -0.02, 0.001, 0, 0),
            new Variable("amp_imag2", -0.07, 0.001, 0, 0), LSKRP, SFKRP, 2);
    Amplitude* Bose_symmetrized_AMP_D = new Amplitude("K*(892)rho(770)_D", new Variable("amp_real3", -0.4, 0.001, 0, 0),
            new Variable("amp_imag3", 0.4, 0.001, 0, 0), LSKRD, SFKRD, 2);

    Amplitude* Bose_symmetrized_KF    = new Amplitude("KF", new Variable("amp_real4", 0.0120787), new Variable("amp_imag4",
            -0.0332525), LSKF, SFKF, 2);
    Amplitude* Bose_symmetrized_KK    = new Amplitude("LSKK", new Variable("amp_real5", 0.0109033),
            new Variable("amp_imag5", -0.00186219), LSKK, SFKK, 2);
    Amplitude* Bose_symmetrized_K1R   = new Amplitude("LSK1R", new Variable("amp_real6", -0.10728),
            new Variable("amp_imag6",  -0.130213), LSK1R, SFK1R, 2);
    Amplitude* Bose_symmetrized_A1R   = new Amplitude("LSA1R", new Variable("amp_real7", 1.0), new Variable("amp_imag7",
            0.0), LSA1R, SFA1R, 2);
    Amplitude* Bose_symmetrized_A1RD  = new Amplitude("LSA1RD", new Variable("amp_real8", -0.94921),
            new Variable("amp_imag8",  -1.73407), LSA1RD, SFA1RD, 2);

    DK3P_DI->amplitudes.push_back(Bose_symmetrized_KF);
    DK3P_DI->amplitudes.push_back(Bose_symmetrized_AMP_S);
    DK3P_DI->amplitudes.push_back(Bose_symmetrized_AMP_P);
    DK3P_DI->amplitudes.push_back(Bose_symmetrized_AMP_D);
    DK3P_DI->amplitudes.push_back(Bose_symmetrized_KK);
    DK3P_DI->amplitudes.push_back(Bose_symmetrized_K1R);
    DK3P_DI->amplitudes.push_back(Bose_symmetrized_A1R);
    DK3P_DI->amplitudes.push_back(Bose_symmetrized_A1RD);

    for(auto res = LSKF.begin(); res != LSKF.end(); ++res) {
        (*res)->setParameterConstantness(true);
    }

    for(auto res = LSKRS.begin(); res != LSKRS.end(); ++res) {
        (*res)->setParameterConstantness(true);
    }

    for(auto res = LSKRP.begin(); res != LSKRP.end(); ++res) {
        (*res)->setParameterConstantness(true);
    }

    for(auto res = LSKRD.begin(); res != LSKRD.end(); ++res) {
        (*res)->setParameterConstantness(true);
    }

    for(auto res = LSKK.begin(); res != LSKK.end(); ++res) {
        (*res)->setParameterConstantness(true);
    }

    for(auto res = LSK1R.begin(); res != LSK1R.end(); ++res) {
        (*res)->setParameterConstantness(true);
    }

    for(auto res = LSA1R.begin(); res != LSA1R.end(); ++res) {
        (*res)->setParameterConstantness(true);
    }

    for(auto res = LSA1RD.begin(); res != LSA1RD.end(); ++res) {
        (*res)->setParameterConstantness(true);
    }

    Variable* constantOne = new Variable("constantOne", 1);
    Variable* constantZero = new Variable("constantZero", 0);

    vector<Variable*> observables;
    vector<Variable*> coefficients;
    vector<Variable*> offsets;

    observables.push_back(m12);
    observables.push_back(m34);
    observables.push_back(cos12);
    observables.push_back(cos34);
    observables.push_back(phi);
    observables.push_back(eventNumber);
    offsets.push_back(constantZero);
    offsets.push_back(constantZero);
    coefficients.push_back(constantOne);

    PolynomialPdf* eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0);
    DPPdf* dp = new DPPdf("test", observables, DK3P_DI, eff, 1e6);

    Variable* constant = new Variable("constant", 0.1);
    Variable* constant2 = new Variable("constant", 1.0);
    vars.clear();
    vars.push_back(constant);
    PolynomialPdf backgr("backgr", m12, vars);
    AddPdf* signal = new AddPdf("signal", constant2, dp, &backgr);

    signal->setData(&currData);
    dp->setDataSize(currData.getNumEvents(), 6);

    FitManager datapdf(signal);
    datapdf.fit();

    return 0;
}
