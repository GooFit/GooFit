
// GooFit stuff
#include <fstream>
#include <goofit/Application.h>
#include <goofit/PDFs/basic/PolynomialPdf.h>
#include <goofit/PDFs/combine/AddPdf.h>
#include <goofit/PDFs/physics/Amp4Body_TD.h>
#include <goofit/PDFs/physics/Lineshapes.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/PDFs/physics/TruthResolution.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>
#include <goofit/FitManager.h>
#include <thrust/count.h>

//#include <goofit/Python.h>
#include <pybind11/embed.h>
namespace py = pybind11;

using namespace std;
using namespace GooFit;

// Constants used in more than one PDF component.
const fptype _mD0       = 1.8645;
const fptype piPlusMass = 0.13957018;
const fptype KmMass     = .493677;

int main(int argc, char **argv) {
    //GooFit::Application app("Time dependent Dalitz plot, 4 particles", argc, argv);


    int trials = 100;
    //app.add_option("-t,--trials", trials, "Number of trials", true);

    //GOOFIT_PARSE(app);

    DecayInfo4t DK3P_DI{Variable("tau", 0.4101, 0.001, 0.3, 0.500),
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
    Variable KstarM{"KstarM", 0.89581, 0.01, 0.8, 1.0};
    Variable KstarW{"KstarW", 0.0474, 0.01, 0.01, 0.2};

    // Variable* f600M  = new Variable("f600M", 0.519, 0.01, 0.75, 0.85);
    // Variable* f600W  = new Variable("f600W", 0.454, 0.01, 0.75, 0.85);
    // Variable* a1M  = new Variable("a1M", 1.23, 0.01, 1.2, 1.3);
    // Variable* a1W  = new Variable("a1W", 0.42, 0.01, 0.37, 0.47);
    // Variable* K1M  = new Variable("K1M", 1.272, 0.01, 1.2, 1.3);
    // Variable* K1W  = new Variable("K1W", 0.09, 0.01, 0.08, 0.1);
    // Variable* K1430M  = new Variable("K1430M", 1.414, 0.01, 1.4, 1.5);
    // Variable* K1430W  = new Variable("K1430W", .29, 0.01, 0.25, 0.35);

    // Spin factors: we have two due to the bose symmetrization of the two pi+
    std::vector<SpinFactor *> SFKRS = {new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S, _mD0, 0, 1, 2, 3),
                                       new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S, _mD0, 3, 1, 2, 0)};

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
        "B_K*(892)rho(770)_S", Variable("amp_real1_B", 1.0), Variable("amp_imag1_B", 0), LSKRS, SFKRS, 2};
    Amplitude Bose_symmetrized_AMP_P_B{
        "B_K*(892)rho(770)_P", Variable("amp_real2_B", -0.145), Variable("amp_imag2_B", 0.86), LSKRP, SFKRP, 2};
    Amplitude Bose_symmetrized_AMP_D_B{
        "B_K*(892)rho(770)_D", Variable("amp_real3_B", 24.343), Variable("amp_imag3_B", 5.329), LSKRD, SFKRD, 2};

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
    Observable eff_weight{"eff_weight",-1000,1000};
    Variable constantOne{"constantOne", 1};
    Variable constantZero{"constantZero", 0};

    //for the master branch don't use eff weights
    vector<Observable> observables{m12, m34, cos12, cos34, phi, eventNumber, dtime, sigmat,eff_weight}; //note order of observables, they are accessed in Amp4Body_TD, in this order!
    //vector<Observable> observables{m12, m34, cos12, cos34, phi, eventNumber, dtime, sigmat}; //note order of observables, they are accessed in Amp4Body_TD, in this order!
    vector<Variable> offsets{constantZero, constantZero};
    vector<Variable> coefficients{constantOne};

    TruthResolution dat;
    PolynomialPdf eff{"constantEff", observables, coefficients, offsets, 0};
   
    
    UnbinnedDataSet fit_data = UnbinnedDataSet(observables);

    double tm12,tm34,tc12,tc34,tphi,tdtime;
    unsigned int fitEvts = 0;
    //read in file
    ifstream mc_file("amp4body_cpp_generated_events.txt");
    if (mc_file.is_open()){
        while(mc_file >> tm12 >> tm34 >> tc12 >> tc34 >> tphi >> tdtime){
            //std::cout << tm12 << " " << tm34 << " " << tc12 << " " << tc34 << " " << tphi << " "<< tdtime << " " << std::endl;

            m12.setValue(tm12);
            m34.setValue(tm34);
            cos12.setValue(tc12);
            cos34.setValue(tc34);
            phi.setValue(tphi);
            dtime.setValue(tdtime);
            sigmat.setValue(0.0);
            eff_weight.setValue(1.0);
            eventNumber.setValue(fitEvts);
            fit_data.addEvent();
            fitEvts++;
        }
    }
    else {
        std::cout << "Cannot open file" << std::endl;
        return 1;
    }
    mc_file.close();

    Amp4Body_TD *dp;
    dp = new Amp4Body_TD("test", observables, DK3P_DI, &dat, &eff, 0, 18000000,true,0);
    printf("current generation seed: %u\n",dp->getGenerationOffset());
    //Amp4Body_TD dp{"test", observables, DK3P_DI, &dat, &eff, 0, 18000000,true};
    std::cout << "Setting special integral flag" << std::endl;
    //dp->set_special_integral(true);
    Variable constant("constant1", 1.0);
    Variable constant2("constant2", 1.0);
    std::vector<Variable> backgrVars = {constant};
    PolynomialPdf backgr("backgr", m12, backgrVars);
    AddPdf signal("signal", constant2, dp, &backgr);
    signal.setData(&fit_data);
    //dp->setData(&fit_data);
    
    dp->setDataSize(fit_data.getNumEvents(),9);
    std::cout << "Fitting data" << std::endl;
    FitManager datapdf(&signal);
    //set special integral
    
    auto fitResults = datapdf.fit();

    std::cout << "calculating pdf values" << std::endl;
    std::vector<std::vector<double>> pdfvals = signal.getCompProbsAtDataPoints();
    std::ofstream pdf_file("amp4body_example_pdf_values.txt");
    for(int i = 0; i < pdfvals[0].size();i++){
        //std::cout <<"total pdf value: "<< pdfvals[0][i]  << "sig pdf value: " << pdfvals[1][i] << " bkg pdf value: " << pdfvals[2][i] << std::endl;
       pdf_file << pdfvals[0][i] << " "  <<  pdfvals[1][i] << " " << pdfvals[2][i] << "\n";
    }
    /*
    //call embedded python here
    py::scoped_interpreter guard{};
    
    auto py_module = py::module_::import("numpy");
    auto hep_ml = py::module_::import("hep_ml");

    auto math = py::module::import("math");
    auto resultobj = math.attr("sqrt")(2);
    double result = resultobj.cast<double>();
    py::print(result);
    //py::print("Hello, World from Python!");
        */
    return 0;
}
