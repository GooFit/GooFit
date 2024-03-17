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
double pi_MASS = 0.139;

double Decay_MASS = B_MASS;
double d1_MASS    = pi_MASS;
double d2_MASS    = pi_MASS;
double d3_MASS    = pi_MASS;

Variable Mother_Mass("DecayParticle_Mass", Decay_MASS);
Variable Daughter1_Mass("DecayProduct_1_Mass", d1_MASS);
Variable Daughter2_Mass("DecayProduct_2_Mass", d2_MASS);
Variable Daughter3_Mass("DecayProduct_3_Mass", d3_MASS);
// For MIPWA
std::vector<fptype> HH_bin_limits; 
std::vector<Variable> pwa_coefs_amp; 
std::vector<Variable> pwa_coefs_phs; 

// Bins for grid normalization
const int bins = 500;

// Dalitz Limits
const fptype s12_min = (d1_MASS + d2_MASS) * (d1_MASS + d2_MASS);
const fptype s12_max = (Decay_MASS - d3_MASS) * (Decay_MASS - d3_MASS);
const fptype s13_min = (d1_MASS + d3_MASS) * (d1_MASS + d3_MASS);
const fptype s13_max = (Decay_MASS - d2_MASS) * (Decay_MASS - d2_MASS);
const fptype s23_min = (d2_MASS + d3_MASS) * (d2_MASS + d3_MASS);
const fptype s23_max = (Decay_MASS - d1_MASS) * (Decay_MASS - d1_MASS);

// Observables
Observable mprime("mprime",0.,1.);
Observable thetaprime("thetaprime",0.,1.);
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


ResonancePdf *loadPWAResonance(std::string name = "pwa_coefs.txt") {
	
        //resonante region < 2GeV
        std::ifstream file(name.c_str());
        double m=0, mag=0, phs=0;
        int index =0;
        while(file >> m >> mag >> phs){
            HH_bin_limits.push_back(m*m);
            Variable va(fmt::format("pwa_coef_{}_real", index), mag ,0.001,0,0);
            Variable vp(fmt::format("pwa_coef_{}_imag", index), phs ,0.001,0,0);
            pwa_coefs_amp.push_back(va);
            pwa_coefs_phs.push_back(vp);
            index++;

        }
        file.close();

        std::cout << "------------------------------------------" << std::endl;
	    std::cout << pwa_coefs_amp.size() << " QMIPWA points loaded!" << std::endl;
	    std::cout << "------------------------------------------" << std::endl;

        for(int i=0; i< HH_bin_limits.size() ; i++){
            std::cout << HH_bin_limits[i]<< " " << pwa_coefs_amp[i].getValue() << " " << pwa_coefs_phs[i].getValue() << std::endl;
        }
				
	    Variable swave_amp_real("swave_real_coef", 1.0);
	    Variable swave_amp_imag("swave_imag_coef", 0.0);



		auto swave = new Resonances::Spline("MIPWA", swave_amp_real, swave_amp_imag, HH_bin_limits, pwa_coefs_amp, pwa_coefs_phs, PAIR_12, true, true);
		
        return swave;
    
} 

Amp3BodySqDP *makesignalpdf(Observable mprime, Observable thetaprime, EventNumber eventNumber, GooPdf *eff = 0) {
    // set up the decay channel
    DecayInfo3 D2KKK;
    D2KKK.motherMass          = B_MASS;
    D2KKK.daug1Mass           = d1_MASS;
    D2KKK.daug2Mass           = d2_MASS;
    D2KKK.daug3Mass           = d3_MASS;
    D2KKK.meson_radius        = 4.0; // GeV^-1
    D2KKK.mother_meson_radius = 5.0; // GeV^-1

    //S-wave
    auto NR = new Resonances::BelleNR("BelleNR",
                                        Variable("BelleNR_amp",0.54*cos(-0.84),0.01,0,0), Variable("BelleNR_img",0.54*sin(-0.84),0.01,0,0), 
                                        Variable("BelleNR_alpha",0.2), PAIR_13, true);

    double f0_980_MASS = 0.965;
    double f0_980_GPP  = 0.2;
    double f0_980_GKK  = 1.0;
    double f0_980_re   = 0.27*cos(-1.59);
    double f0_980_img  = 0.27*sin(-1.59);

    Variable v_f0_980_Mass("f0_980_MASS", f0_980_MASS);
    Variable v_f0_980_GPP("f0_980_GPP", f0_980_GPP);
    Variable v_f0_980_GKK("f0_980_GKK", f0_980_GKK);
    Variable v_f0_980_real("f0_980_REAL", f0_980_re, 0.01, 0, 0);
    Variable v_f0_980_img("f0_980_IMAG", f0_980_img, 0.01, 0, 0);

    auto f0_980 = new Resonances::FLATTE(
        "f0_980", v_f0_980_real, v_f0_980_img, v_f0_980_Mass, v_f0_980_GPP, v_f0_980_GKK, 0, PAIR_13 ,  true);
    //P-wave

    double rho770_MASS  = 0.77526;
    double rho770_WIDTH = 0.1478;
    double rho770_amp   = 1.;
    double rho770_img   = 0.;

    Variable v_rho770_Mass("rho770_MASS", rho770_MASS);
    Variable v_rho770_Width("rho770_WIDTH", rho770_WIDTH);
    Variable v_rho770_real("rho770_REAL", rho770_amp, 0.01, 0, 0);
    Variable v_rho770_img("rho770_IMAG", rho770_img, 0.01, 0, 0);

    auto rho770 = new Resonances::RBW(
        "rho770", v_rho770_real, v_rho770_img, v_rho770_Mass, v_rho770_Width, 1, PAIR_13, true, true);

    double rho1450_MASS  = 1.465;
    double rho1450_WIDTH = 0.4;
    double rho1450_amp   = 0.37*cos(1.99);
    double rho1450_img   = 0.37*sin(1.99);

    Variable v_rho1450_Mass("rho1450_MASS", rho1450_MASS);
    Variable v_rho1450_Width("rho1450_WIDTH", rho1450_WIDTH);
    Variable v_rho1450_real("rho1450_REAL", rho1450_amp,0.01,0,0);
    Variable v_rho1450_img("rho1450_IMAG", rho1450_img,0.01,0,0);


    auto rho1450 = new Resonances::RBW(
        "rho1450", v_rho1450_real, v_rho1450_img, v_rho1450_Mass, v_rho1450_Width, 1, PAIR_13, true, true);

    // D-wave
    double f2_1270_MASS  = 1.2751;
    double f2_1270_WIDTH = 0.1851;
    double f2_1270_amp   = 1.0;//0.53*cos(1.39);
    double f2_1270_img   = 0.0;//0.53*sin(1.39);

    Variable v_f2_1270_Mass("f2_1270_MASS", f2_1270_MASS);
    Variable v_f2_1270_Width("f2_1270_WIDTH", f2_1270_WIDTH);
    Variable v_f2_1270_real("f2_1270_REAL", f2_1270_amp);//,0.001,0,0);
    Variable v_f2_1270_img("f2_1270_IMAG", f2_1270_img);//,0.001,0,0);

    auto f2_1270 = new Resonances::RBW(
        "f2", v_f2_1270_real, v_f2_1270_img, v_f2_1270_Mass, v_f2_1270_Width, 2, PAIR_13, true, true);

    //swave
    auto swave = loadPWAResonance("PWAFile.bin");

    
    std::vector<ResonancePdf *> vec_resonances;
    vec_resonances.push_back(rho770);
    vec_resonances.push_back(rho1450);
    vec_resonances.push_back(f0_980);
    vec_resonances.push_back(f2_1270);
    vec_resonances.push_back(NR);

    //vec_resonances.push_back(swave);

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
    auto flags_val(0.);

    printf("NEntries = %f \n",tree->GetEntries());

    tree->SetBranchAddress("mPrime", &mprime_val);
    tree->SetBranchAddress("thPrime", &thetaprime_val);

   
    size_t j = 0;
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(16); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.,1.);
    for(size_t i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);

        mprime.setValue(mprime_val);
        if(dis(gen)>0.5)
            thetaprime.setValue(thetaprime_val);
        else
            thetaprime.setValue(1.-thetaprime_val);
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

    double _mprime, _thetaprime,_s12,_s13,_s23,_jac;
    auto f     = new TFile(name.c_str(), "recreate");
    auto t     = new TTree("genResults", "");
    auto b_mprime = t->Branch("mPrime", &_mprime, "mPrime/D");
    auto b_thetaprime = t->Branch("thPrime", &_thetaprime, "thPrime/D");
    t->Branch("s12", &_s12, "s12/D");
    t->Branch("s13", &_s13, "s13/D");
    t->Branch("s23", &_s23, "s23/D");
    t->Branch("jac", &_jac, "jac/D");
    

    for(int i = 0; i < toyMC.getNumEvents(); i++) {
        toyMC.loadEvent(i);
        t->GetEntry(i);
        _mprime = mprime.getValue();
        _thetaprime = thetaprime.getValue();
        // if(_thetaprime>0.5)
        //     _thetaprime = 1.0 -_thetaprime;
        _s12 = pow(calc_m12(_mprime, Decay_MASS, d1_MASS, d2_MASS, d3_MASS),2);
        _s13 = pow(calc_m13(sqrt(_s12), cos(_thetaprime*M_PI), Decay_MASS, d1_MASS, d2_MASS, d3_MASS),2);
        _s23 = Decay_MASS*Decay_MASS + d1_MASS*d1_MASS + d2_MASS*d2_MASS + d3_MASS*d3_MASS - _s13 - _s12;
        _jac = calc_SqDp_Jacobian(_mprime, _thetaprime, Decay_MASS, d1_MASS, d2_MASS, d3_MASS);
        t->Fill();
    }
    t->Write("", TObject::kOverwrite);
    f->Write("", TObject::kOverwrite);
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
    //datapdf.printParams();
    output = fmt::format("Fit/{0}/fit_result_fitted.txt", name.c_str());
    writeToFile(totalPdf, output.c_str());

    return signal;
}

int main(int argc, char **argv) {
    GooFit::Application app{"Genfit", argc, argv};

    std::string input_data_name = "gen-3K.root";
    std::string fit_name        = "Fit";
    std::string acc_file        = "acc_hist_0_Smoothed.root";
    std::string bkg_file        = "bkg_hist_0_Smoothed.root";
    std::string toyName         = "MC.root";
    bool save_toy               = false;
    bool is_toy                 = false;
    bool no_acc_and_bkg         = true;
    size_t Nevents              = 100000;

    auto fit = app.add_subcommand("fit", "fit data");
    fit->add_option("-f,--file", input_data_name, "name_of_file.root");
    fit->add_option("-t,--isToy", is_toy, "Get toyData for fit");
    fit->add_option("-s,--saveToy", save_toy, "save toy in root file");
    fit->add_option("-n,--fitName", fit_name, "name of this fit(useful to save results)");
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

    auto makeToyMCB = app.add_subcommand("makeToyMCB", "make toy using MC Booster");
    makeToyMCB->add_option("-e,--nevents", Nevents, "number of events");
    makeToyMCB->add_option("-n,--name", toyName, "output_toy_name.root");
    makeToyMCB->add_option("-s,--saveToy", save_toy, "save toy in root file");
    makeToyMCB->add_option("-a,--acc", acc_file, "name of acc file");
    makeToyMCB->add_option("-b,--bkg", bkg_file, "name of bkg file");
    makeToyMCB->add_option("-d,--disable-acc-bkg", no_acc_and_bkg, "disable-acc-bkg");

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
    GooPdf *totalpdf = nullptr;

    auto signal = makesignalpdf(mprime, thetaprime, eventNumber, efficiency);
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
        auto  fullName = fmt::format("Fit/{0}/{1}",fit_name, "output_plots.root");
        SqDalitzPlotter dplotter{totalpdf, signal};
        dplotter.Plot(fullName, &data);
    
        data.clear();
        dplotter.fillDataSetMC(data, Nevents);
        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << data.getNumEvents() << " events was generated!" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
        fullName = fmt::format("Fit/{0}/{1}",fit_name, toyName);
        
        //to_root(data, fullName);
        std::cout << toyName << " root file was saved in MC folder" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
        //std::cout << "Norm: " << signal->normalize() << std::endl;
        auto frac = signal->fit_fractions(true);

        std::vector<fptype> pwa_coefs_amp_fitted; 
        std::vector<fptype> pwa_coefs_phs_fitted; 

        for(int i=0; i<HH_bin_limits.size() ; i++){
            auto real = pwa_coefs_amp[i].getValue();
            auto imag = pwa_coefs_phs[i].getValue();
            auto mag = sqrt(real*real + imag*imag);
            auto phs = atan2(imag,real);
            pwa_coefs_amp_fitted.push_back(mag);
            pwa_coefs_phs_fitted.push_back(phs);
        }

        TGraph gr_mag(HH_bin_limits.size(),HH_bin_limits.data(),pwa_coefs_amp_fitted.data());
        TGraph gr_phs(HH_bin_limits.size(),HH_bin_limits.data(),pwa_coefs_phs_fitted.data());

        TCanvas c;
        gr_mag.Draw("APL*");
        c.SaveAs("pwa_coefs_amp_fitted.png");
        gr_phs.Draw("APL*");
        c.SaveAs("pwa_coefs_phs_fitted.png");
       
    }

    if(*makeToyMCB) {
        
        
            auto name= fmt::format("MC/{0}",toyName);
            double _s12, _s13,_s23, mp, th, _jac;
            double fg, wt;
            auto f = new TFile(name.c_str(),"recreate");
            auto t = new TTree("genResults","toyMC");
            auto b_mp  = t->Branch("mPrime",&mp,"mPrime/D");
            auto b_th  = t->Branch("thPrime",&th,"thPrime/D");	
            auto b_s12 = t->Branch("s12",&_s12,"s12/D");
            auto b_s13 = t->Branch("s13",&_s13,"s13/D");
            auto b_s23 = t->Branch("s23",&_s23,"s23/D");
            auto b_flags = t->Branch("flags",&fg,"flags/D");
            auto b_wt = t->Branch("weights",&wt,"weights/D");
            auto b_jac = t->Branch("jac",&_jac,"jac/D");
            
            std::srand(std::time(nullptr)); // use current time as seed for random generator
           
            int accepted = 0;
            while(t->GetEntries() < Nevents){
            
                int random_value = std::rand();
                auto tuple = signal->GenerateSig(1000000,random_value);
                auto variables = std::get<1>(tuple);
                auto weights   = std::get<2>(tuple);
                auto flags     = std::get<3>(tuple);
                 
                std::cout << "events = " << t->GetEntries() << "\n";
                for(int i = 0; i < weights.size(); i++){
                        if(flags[i]==1){
                            t->GetEntry(accepted);
                            _s12 = (*(variables[4]))[i];
                            _s13 = (*(variables[2]))[i];
                            _s23 = (*(variables[3]))[i];
                            mp = (*(variables[0]))[i];
                            th = (*(variables[1]))[i];
                            fg = flags[i];
                            wt = weights[i];
                            _jac = calc_SqDp_Jacobian(mp, th, Decay_MASS, d1_MASS, d2_MASS, d3_MASS);
                            t->Fill();   
                            accepted++;
                        }
                }
            }

            t->Write("",TObject::kOverwrite);
            f->Write();
            f->Close();
            std::cout << "------------------------------------------" << std::endl;
            std::cout << "toyMC --> " << name.c_str() << " was saved!" << std::endl;
            std::cout << "------------------------------------------" << std::endl;

            auto frac = signal->fit_fractions(true);

     

            //Plot S-wave
            
        
    }




}
