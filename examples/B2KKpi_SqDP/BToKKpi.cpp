#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TTree.h>
#include <TROOT.h>
#include <TMinuit.h>
#include <TNtuple.h>
#include <TComplex.h>
#include <TFile.h>
#include <TStyle.h>
#include <TH2Poly.h>
#include <TGraphErrors.h>

// Minuit
#include <Minuit2/MnStrategy.h>
#include <Minuit2/Minuit2Minimizer.h>

// System stuff
#include <CLI/Timer.hpp>
#include <fstream>
#include <string>
#include <time.h>

// GooFit stuff
#include <goofit/Application.h>
#include <goofit/BinnedDataSet.h>
#include <goofit/FitManager.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/basic/SmoothHistogramPdf.h>
#include <goofit/PDFs/basic/PolynomialPdf.h>
#include <goofit/PDFs/combine/AddPdf.h>
#include <goofit/PDFs/combine/ProdPdf.h>
#include <goofit/PDFs/physics/Amp3BodySqDP.h>
#include <goofit/PDFs/physics/DalitzVetoPdf.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>
#include <goofit/PDFs/physics/SqDalitzPlotter.h>
#include <goofit/FunctorWriter.h>

// Matrix
#include <Eigen/Core>
#include <Eigen/LU>

#include <thrust/count.h>


#define torad(x) (x * M_PI / 180.)

using namespace std;
using namespace GooFit;
using namespace ROOT;

// B_MASS (PDG)
double B_MASS = 5.27934;
double k_MASS = 0.493677;
double pi_MASS = 0.13957039;

double Decay_MASS = B_MASS;
double d1_MASS    = k_MASS;
double d2_MASS    = pi_MASS;
double d3_MASS    = k_MASS;

Variable Mother_Mass("DecayParticle_Mass", Decay_MASS);
Variable Daughter1_Mass("DecayProduct_1_Mass", d1_MASS);
Variable Daughter2_Mass("DecayProduct_2_Mass", d2_MASS);
Variable Daughter3_Mass("DecayProduct_3_Mass", d3_MASS);

// Bins for grid normalization
const int bins = 1000;

// Dalitz Limits
const fptype s12_min = (d1_MASS + d2_MASS) * (d1_MASS + d2_MASS);
const fptype s12_max = (Decay_MASS - d3_MASS) * (Decay_MASS - d3_MASS);
const fptype s13_min = (d1_MASS + d3_MASS) * (d1_MASS + d3_MASS);
const fptype s13_max = (Decay_MASS - d2_MASS) * (Decay_MASS - d2_MASS);
const fptype s23_min = (d2_MASS + d3_MASS) * (d2_MASS + d3_MASS);
const fptype s23_max = (Decay_MASS - d1_MASS) * (Decay_MASS - d1_MASS);

// Observables
Observable mprime("mprime",0,1);
Observable thetaprime("thetaprime",0,1);
Observable s23("s23", s23_min, s23_max);
EventNumber eventNumber("eventNumber");

GooPdf *polyEff(Observable s12, Observable s13) {
    vector<Variable> offsets;
    vector<Observable> observables;
    vector<Variable> coefficients;
    Variable constantOne("c1", 1);
    Variable constantZero("c0", 0);
    observables.push_back(mprime);
    observables.push_back(thetaprime);
    offsets.push_back(constantZero);
    offsets.push_back(constantZero);
    coefficients.push_back(constantOne);
    PolynomialPdf *eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0); // No efficiency

    return eff;
}

GooPdf *Vetos(std::string name) {
    // for Dzero veto
    fptype DMin31 = 1.83*1.83;
	fptype DMax31 = 2.00*2.00;

    VetoInfo veto13(Variable("veto13_min", DMin31), Variable("veto13_max", DMax31), PAIR_13);
    VetoInfo veto23(Variable("veto23_min", DMin31), Variable("veto23_max", DMax31), PAIR_23);

    vector<VetoInfo> vetos;
    vetos.push_back(veto13);
    vetos.push_back(veto23);


    return new DalitzVetoPdf(name, mprime, thetaprime, Mother_Mass, Daughter1_Mass, Daughter2_Mass, Daughter3_Mass, vetos);
}

Amp3BodySqDP *makesignalpdf(Observable mprime, Observable thetaprime, EventNumber eventNumber, GooPdf *eff = 0) {
    // set up the decay channel
    DecayInfo3 D2KKK;
    D2KKK.motherMass          = B_MASS;
    D2KKK.daug1Mass           = d1_MASS;
    D2KKK.daug2Mass           = d2_MASS;
    D2KKK.daug3Mass           = d3_MASS;
    D2KKK.meson_radius        = 4.0; // GeV^-1
    D2KKK.mother_meson_radius = 4.0; // GeV^-1
    //PAIR_12 = K+pi+
    //PAIR_13 = K+K-
    //PAIR_23 = pi+K-
    //Laura (1 => m_23, 2 => m_13, 3 => m_12), i.e. the bachelor track number.
    //Kst(1430)
    double Ks1430_MASS  = 1.425;
    double Ks1430_WIDTH = 0.270;
    double Ks1430_amp   = 1.;
    double Ks1430_img   = 0.;

    Variable v_Ks1430_Mass("Ks1430_MASS", Ks1430_MASS);
    Variable v_Ks1430_Width("Ks1430_WIDTH", Ks1430_WIDTH);
    Variable v_Ks1430_real("Ks1430_REAL", Ks1430_amp,0.01,0,0);
    Variable v_Ks1430_img("Ks1430_IMAG", Ks1430_img,0.01,0,0);

    auto Ks1430 = new Resonances::RBW(
        "Ks1430", v_Ks1430_real, v_Ks1430_img, v_Ks1430_Mass, v_Ks1430_Width, 0, PAIR_23, false, false);

     //chi_c0
    double Chic0_MASS  = 3.41475;
    double Chic0_WIDTH = 0.0105;
    double Chic0_amp   = 1.;
    double Chic0_img   = 0.;

    Variable v_Chic0_Mass("Chic0_MASS", Chic0_MASS);
    Variable v_Chic0_Width("Chic0_WIDTH", Chic0_WIDTH);
    Variable v_Chic0_real("Chic0_REAL", Chic0_amp,0.01,0,0);
    Variable v_Chic0_img("Chic0_IMAG", Chic0_img,0.01,0,0);

    auto Chic0 = new Resonances::RBW(
        "Chic0", v_Chic0_real, v_Chic0_img, v_Chic0_Mass, v_Chic0_Width, 0, PAIR_13, false, false);

    double Ks892_MASS  = 0.89581;
    double Ks892_WIDTH = 0.0474;
    double Ks892_amp   = 1.;
    double Ks892_img   = 0.;

    //K*(892)
    Variable v_Ks892_Mass("Ks892_MASS", Ks892_MASS);
    Variable v_Ks892_Width("Ks892_WIDTH", Ks892_WIDTH);
    Variable v_Ks892_real("Ks892_REAL", Ks892_amp);
    Variable v_Ks892_img("Ks892_IMAG", Ks892_img);

    auto Ks892 = new Resonances::RBW(
        "Ks892", v_Ks892_real, v_Ks892_img, v_Ks892_Mass, v_Ks892_Width, 1, PAIR_23, false, false);

    //phi(1020)
    double phi1020_MASS  = 1.019461;
    double phi1020_WIDTH = 0.004266;
    double phi1020_amp   = 1.;
    double phi1020_img   = 0.;

    Variable v_phi1020_Mass("phi1020_MASS", phi1020_MASS);
    Variable v_phi1020_Width("phi1020_WIDTH", phi1020_WIDTH);
    Variable v_phi1020_real("phi1020_REAL", phi1020_amp,0.01,0,0);
    Variable v_phi1020_img("phi1020_IMAG", phi1020_img,0.01,0,0);

    auto phi1020 = new Resonances::RBW(
        "phi1020", v_phi1020_real, v_phi1020_img, v_phi1020_Mass, v_phi1020_Width, 1, PAIR_13, false, false);
    
    //rho(1450)
    double rho1450_MASS  = 1.465;
    double rho1450_WIDTH = 0.4;
    double rho1450_amp   = 1.;
    double rho1450_img   = 0.;

    Variable v_rho1450_Mass("rho1450_MASS", rho1450_MASS);
    Variable v_rho1450_Width("rho1450_WIDTH", rho1450_WIDTH);
    Variable v_rho1450_real("rho1450_REAL", rho1450_amp,0.01,0,0);
    Variable v_rho1450_img("rho1450_IMAG", rho1450_img,0.01,0,0);

    auto rho1450 = new Resonances::RBW(
        "rho1450", v_rho1450_real, v_rho1450_img, v_rho1450_Mass, v_rho1450_Width, 1, PAIR_13, false, false);

    //rho(1700)
    double rho1700_MASS  = 1.720;
    double rho1700_WIDTH = 0.250;
    double rho1700_amp   = 1.;
    double rho1700_img   = 0.;

    Variable v_rho1700_Mass("rho1700_MASS", rho1700_MASS);
    Variable v_rho1700_Width("rho1700_WIDTH", rho1700_WIDTH);
    Variable v_rho1700_real("rho1700_REAL", rho1700_amp,0.01,0,0);
    Variable v_rho1700_img("rho1700_IMAG", rho1700_img,0.01,0,0);

    auto rho1700 = new Resonances::RBW(
        "rho1700", v_rho1700_real, v_rho1700_img, v_rho1700_Mass, v_rho1700_Width, 1, PAIR_13, false, false);

    //f2(1270)
    double f2_1270_MASS  = 1.2751;
    double f2_1270_WIDTH = 0.1851;
    double f2_1270_amp   = 0.;
    double f2_1270_img   = 1.;

    Variable v_f2_1270_Mass("f2_1270_MASS", f2_1270_MASS);
    Variable v_f2_1270_Width("f2_1270_WIDTH", f2_1270_WIDTH);
    Variable v_f2_1270_real("f2_1270_REAL", f2_1270_amp, 0.01, 0, 0);
    Variable v_f2_1270_img("f2_1270_IMAG", f2_1270_img, 0.01, 0, 0);

    auto f2_1270 = new Resonances::RBW(
        "f2_1270", v_f2_1270_real, v_f2_1270_img, v_f2_1270_Mass, v_f2_1270_Width, 2, PAIR_13, false, false);
    
    //f'2(1525)
    double f2p_1525_MASS  = 1.525;
    double f2p_1525_WIDTH = 0.073;
    double f2p_1525_amp   = 0.;
    double f2p_1525_img   = 1.;

    Variable v_f2p_1525_Mass("f2p_1525_MASS", f2p_1525_MASS);
    Variable v_f2p_1525_Width("f2p_1525_WIDTH", f2p_1525_WIDTH);
    Variable v_f2p_1525_real("f2p_1525_REAL", f2p_1525_amp, 0.01, 0, 0);
    Variable v_f2p_1525_img("f2p_1525_IMAG", f2p_1525_img, 0.01, 0, 0);

    auto f2p_1525 = new Resonances::RBW(
        "f2p_1525", v_f2p_1525_real, v_f2p_1525_img, v_f2p_1525_Mass, v_f2p_1525_Width, 2, PAIR_13, false, false);


    auto nonres = new Resonances::NonRes("NonRes",Variable("NR_REAL",1.,0.01,0,0),Variable("NR_IMAG",0.,0.01,0,0));

    auto PolarFFNR = new Resonances::PolarFFNR("PolarFFNR",Variable("FFNR_REAL",1.,0.01,0,0),Variable("FFNR_IMAG",0.,0.01,0,0),Variable("lambda",1.18),PAIR_23,false);

    std::vector<Variable> ScatteringCoefs;
    ScatteringCoefs.push_back(Variable("B1Val",23.6));
    ScatteringCoefs.push_back(Variable("B2Val",29.4));
    ScatteringCoefs.push_back(Variable("B3Val",0.6));
    ScatteringCoefs.push_back(Variable("C1Val",34.39));
    ScatteringCoefs.push_back(Variable("C2Val",4.4));
    ScatteringCoefs.push_back(Variable("C3Val",-32.9));
    ScatteringCoefs.push_back(Variable("C4Val",-16.));
    ScatteringCoefs.push_back(Variable("C5Val",7.4));
    ScatteringCoefs.push_back(Variable("D0Val",0.59));
    ScatteringCoefs.push_back(Variable("D1Val",-0.38));
    ScatteringCoefs.push_back(Variable("D2Val",0.12));
    ScatteringCoefs.push_back(Variable("D3Val",-0.09));
    ScatteringCoefs.push_back(Variable("F1Val",-0.043));
    ScatteringCoefs.push_back(Variable("F2Val",-0.008));
    ScatteringCoefs.push_back(Variable("F3Val",-0.28));
    ScatteringCoefs.push_back(Variable("F4Val",0.026));

    auto Scattering = new Resonances::Rescattering2("Scattering",Variable("Scat_REAL",1.,0.01,0,0),Variable("Scat_IMAG",0.,0.01,0,0),ScatteringCoefs,PAIR_13,false);


    std::vector<ResonancePdf *> vec_resonances;

    vec_resonances.push_back(Ks892);
    vec_resonances.push_back(phi1020);
    vec_resonances.push_back(f2_1270);
    vec_resonances.push_back(Ks1430);
    vec_resonances.push_back(rho1450);
    vec_resonances.push_back(f2p_1525);
    vec_resonances.push_back(rho1700);
    vec_resonances.push_back(Chic0);
    vec_resonances.push_back(PolarFFNR);
    vec_resonances.push_back(Scattering);

    D2KKK.resonances = vec_resonances;

    return new Amp3BodySqDP("signalPDF", mprime, thetaprime, eventNumber, D2KKK, eff);
}

void getData(std::string toyFileName, GooFit::Application &app, DataSet &data, bool toy) {
    // load data in a GooFit::dataset

    auto obs               = data.getObservables();
    Observable mprime         = obs.at(0);
    Observable thetaprime         = obs.at(1);
    Observable eventNumber = obs.at(2);

    auto openRoot = new TFile(toyFileName.c_str());
    auto tree     = (TTree *)openRoot->Get("genResults");
    auto mprime_val(0.);
    auto thetaprime_val(0.);
    int charge(0);

    printf("NEntries = %f \n",tree->GetEntries());

    tree->SetBranchAddress("mPrime", &mprime_val);
    tree->SetBranchAddress("thPrime", &thetaprime_val);
    tree->SetBranchAddress("charge", &charge);
   
    size_t j = 0;
    for(size_t i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        if(charge<0)
            continue;

        mprime.setValue(mprime_val);
        thetaprime.setValue(thetaprime_val);
        eventNumber.setValue(data.getNumEvents());
        data.addEvent();
        if(j < 10)
            printf("[%d] = (%f , %f)\n", i, mprime.getValue(), thetaprime.getValue());
        j++;
      
    }
}

void to_root(UnbinnedDataSet &toyMC, std::string name) {
    // save GooFit::Dataset in a root file
    auto obs               = toyMC.getObservables();
    Observable mprime         = obs.at(0);
    Observable thetaprime         = obs.at(1);
    Observable eventNumber = obs.at(2);

    double _mprime, _thetaprime,_s12,_s13,_s23;
    auto f     = new TFile(name.c_str(), "recreate");
    auto t     = new TTree("DecayTree", "");
    auto b_mprime = t->Branch("mprime", &_mprime, "mprime/D");
    auto b_thetaprime = t->Branch("thetaprime", &_thetaprime, "thetaprime/D");
    t->Branch("s12", &_s12, "s12/D");
    t->Branch("s13", &_s13, "s13/D");
    t->Branch("s23", &_s23, "s23/D");
    

    for(int i = 0; i < toyMC.getNumEvents(); i++) {
        toyMC.loadEvent(i);
        t->GetEntry(i);
        _mprime = mprime.getValue();
        _thetaprime = thetaprime.getValue();
        _s12 = pow(calc_m12(_mprime, Decay_MASS, d1_MASS, d2_MASS, d3_MASS),2);
        _s13 = pow(calc_m13(sqrt(_s12), cos(_thetaprime*M_PI), Decay_MASS, d1_MASS, d2_MASS, d3_MASS),2);
        _s23 = Decay_MASS*Decay_MASS + d1_MASS*d1_MASS + d2_MASS*d2_MASS + d3_MASS*d3_MASS - _s13 - _s12;
        t->Fill();
    }
    t->Write("", TObject::kOverwrite);
    f->Write();
    f->Close();
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "toyMC --> " << name.c_str() << " was saved!" << std::endl;
    std::cout << "------------------------------------------" << std::endl;
}

Amp3BodySqDP *runFit(GooPdf *totalPdf, Amp3BodySqDP *signal, UnbinnedDataSet *data, std::string name) {
    // This function run the data fit

    // Setting data and EventNumber
    totalPdf->setData(data);
    signal->setDataSize(data->getNumEvents());

    auto output = fmt::format("Fit/{0}/fit_result_before.txt", name.c_str());
    writeToFile(totalPdf, output.c_str());

    // Fitter (it uses ROOT::FunctionMinimum API)
    FitManager datapdf(totalPdf);
    datapdf.setVerbosity(2);
    datapdf.setMaxCalls(200000);

    // Start fit
    auto func_min = datapdf.fit();
    datapdf.printParams();
    output = fmt::format("Fit/{0}/fit_result_fitted.txt", name.c_str());
    writeToFile(totalPdf, output.c_str());

    return signal;
}

int main(int argc, char **argv) {
    GooFit::Application app{"Genfit", argc, argv};

    std::string input_data_name = "input.root";
    std::string fit_name        = "Fit";
    std::string acc_file        = "acc_hist_0_Smoothed.root";
    std::string bkg_file        = "bkg_hist_0_Smoothed.root";
    std::string toyName         = "MC.root";
    bool save_toy               = false;
    bool is_toy                 = false;
    bool no_acc_and_bkg         = true;
    size_t Nevents              = 100875;

    auto fit = app.add_subcommand("fit", "fit data");
    fit->add_option("-f,--file", input_data_name, "name_of_file.root");
    fit->add_option("-t,--isToy", is_toy, "Get toyData for fit");
    fit->add_option("-s,--saveToy", save_toy, "save toy in root file");
    fit->add_option("-n,--fitName", fit_name, "name of this fit(useful to save results)")->required(true);
    fit->add_option("-a,--acc", acc_file, "name of acc file");
    fit->add_option("-b,--bkg", bkg_file, "name of bkg file");
    fit->add_option("-d,--disable-acc-bkg", no_acc_and_bkg, "disable-acc-bkg");

    auto makeToy = app.add_subcommand("makeToy", "make a toy");
    makeToy->add_option("-e,--nevents", Nevents, "number of events");
    makeToy->add_option("-n,--name", toyName, "output_toy_name.root");
    makeToy->add_option("-s,--saveToy", save_toy, "save toy in root file");
    makeToy->add_option("-a,--acc", acc_file, "name of acc file");
    makeToy->add_option("-b,--bkg", bkg_file, "name of bkg file");
    makeToy->add_option("-d,--disable-acc-bkg", no_acc_and_bkg, "disable-acc-bkg");

    app.require_subcommand();

    GOOFIT_PARSE(app);

    // Make the MC directory if it does not exist
    std::string command = "mkdir -p MC";
    if(system(command.c_str()) != 0)
        throw GooFit::GeneralError("Making `MC` directory failed");

    // Make the Fit directory if it does not exist
    command = "mkdir -p Fit";
    if(system(command.c_str()) != 0)
        throw GooFit::GeneralError("Making `Fit` directory failed");

    mprime.setNumBins(bins);
    thetaprime.setNumBins(bins);

    auto efficiency  = polyEff(mprime, thetaprime);
    auto veto = Vetos("veto");
    auto vetoEff = new ProdPdf("vetoEff",{efficiency,veto}) ;
    GooPdf *totalpdf = nullptr;

    auto signal = makesignalpdf(mprime, thetaprime, eventNumber, vetoEff);
    totalpdf    = new ProdPdf("totalpdf", {signal});

    


    if(*makeToy) {
        SqDalitzPlotter dplotter{totalpdf, signal};
        UnbinnedDataSet data({mprime, thetaprime, eventNumber});
        dplotter.fillDataSetMC(data, Nevents);
        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << data.getNumEvents() << " events was generated!" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
        if(save_toy) {
            auto fullName = fmt::format("MC/{0}", toyName);
            to_root(data, fullName);
            std::cout << toyName << " root file was saved in MC folder" << std::endl;
            std::cout << "----------------------------------------------------------" << std::endl;
        }

        
        signal->setDataSize(data.getNumEvents());
        totalpdf->setData(&data);
        
        auto frac = signal->fit_fractions(true);
       
    

        return 0;
    }

    // if(*makeToy) {
    //     long ngenerated = 0;
    //     uint generationOffset = 0;
    //     uint batchSize            = 2000000;
    //     auto name= fmt::format("MC/{0}",toyName);
    //     auto f = new TFile(name.c_str(),"recreate");
    //     auto t = new TTree("DecayTree","toyMC");
    //     double _s12, _s13,_s23;
    //     auto b_s12 = t->Branch("mprime",&_s12,"s12/D");
    //     auto b_s13 = t->Branch("thetaprime",&_s13,"s13/D");
    //     auto b_s23 = t->Branch("s23",&_s23,"s23/D");

    //     while(ngenerated<Nevents){

            
    //         printf("ngenerated = %d \n",ngenerated);

    //         signal->setGenerationOffset(generationOffset);
    //         auto tuple = signal->GenerateSig(batchSize,721863);
       
    //         auto variables = std::get<1>(tuple);
    //         auto flags     = std::get<3>(tuple);
    //         long accepted  = thrust::count_if(flags.begin(), flags.end(), thrust::identity<bool>());
            
    //         for(int i = 0; i < flags.size(); i++){
    //             if(ngenerated<Nevents && flags[i] == 1) {
    //                 t->GetEntry(i+ngenerated);
    //                 _s12 = (*(variables[0]))[i];
    //                 _s13 = (*(variables[2]))[i];
    //                 _s23 = (*(variables[1]))[i];
                   
    //                 t->Fill();
    //             }
    //         }

    //         generationOffset += batchSize;
    //         ngenerated+=accepted;
    //         delete variables[0];
    //         delete variables[1];
    //         delete variables[2];
        
    //     }

    //     t->Write("",TObject::kOverwrite);
    //     f->Write("",TObject::kOverwrite);
    //     f->Close();
    //     std::cout << "------------------------------------------" << std::endl;
    //     std::cout << "toyMC --> " << name.c_str() << " was saved!" << std::endl;
    //     std::cout << "------------------------------------------" << std::endl;
    //     return 0;
    // }

    if(*fit) {
        auto command = fmt::format("mkdir -p Fit/{0}", fit_name);
        if(system(command.c_str()) != 0)
            throw GooFit::GeneralError("Making directory failed");

        std::cout << "------------------------------------------" << std::endl;
        std::cout << "Reading file --> " << input_data_name << std::endl;
        std::cout << "------------------------------------------" << std::endl;
        UnbinnedDataSet data({mprime, thetaprime, eventNumber});
        getData(input_data_name, app, data, is_toy);
        std::cout << "------------------------------------------" << std::endl;
        std::cout << "Num Entries Loaded =  " << data.getNumEvents() << std::endl;
        std::cout << "------------------------------------------" << std::endl;
        auto output_signal = runFit(totalpdf, signal, &data, fit_name);
        data.clear();
        SqDalitzPlotter dplotter{totalpdf, signal};
        dplotter.fillDataSetMC(data, Nevents);
        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << data.getNumEvents() << " events was generated!" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
        auto fullName = fmt::format("MC/{0}", toyName);
        to_root(data, fullName);
        std::cout << toyName << " root file was saved in MC folder" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;

   
        signal->setDataSize(data.getNumEvents());
        totalpdf->setData(&data);
    
        auto frac = signal->fit_fractions(true);
       
    }
}