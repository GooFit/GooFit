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
#include <goofit/PDFs/basic/BinTransformPdf.h>
#include <goofit/PDFs/combine/MappedPdf.h>
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
#include <Eigen/Dense>

using namespace Eigen;

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
Observable category("category", 0., 1.);
EventNumber eventNumber("eventNumber");

GooPdf *polyEff(Observable s12, Observable s13,  bool charge_plus=true ) {

    std::string Charge = charge_plus ? "_Plus" : "_Minus";
    vector<Variable> offsets;
    vector<Observable> observables;
    vector<Variable> coefficients;
    Variable constantOne("c1"+Charge, 1);
    Variable constantZero("c0"+Charge, 0);
    observables.push_back(mprime);
    observables.push_back(thetaprime);
    offsets.push_back(constantZero);
    offsets.push_back(constantZero);
    coefficients.push_back(constantOne);
    PolynomialPdf *eff = new PolynomialPdf("constantEff"+Charge, observables, coefficients, offsets, 0); // No efficiency

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

Amp3BodySqDP *makesignalpdf(Observable mprime, Observable thetaprime, EventNumber eventNumber, GooPdf *eff = 0, bool charge_plus=true) {
    // set up the decay channel
    DecayInfo3 D2KKK;
    D2KKK.motherMass          = B_MASS;
    D2KKK.daug1Mass           = d1_MASS;
    D2KKK.daug2Mass           = d2_MASS;
    D2KKK.daug3Mass           = d3_MASS;
    D2KKK.meson_radius        = 4.0; // GeV^-1
    D2KKK.mother_meson_radius = 5.0; // GeV^-1

    std::string Charge = charge_plus ? "_Plus" : "_Minus";



    //S-wave
    auto NR = new Resonances::BelleNR("BelleNR"+Charge,
                                        Variable("BelleNR_amp"+Charge,0.36042993,0.01,0,0), 
                                        Variable("BelleNR_img"+Charge,-0.40210728,0.01,0,0), 
                                        Variable("BelleNR_alpha"+Charge,0.2), PAIR_13, true);

    double f0_980_MASS = 0.965;
    double f0_980_GPP  = 0.2;
    double f0_980_GKK  = 1.0;
    double f0_980_re   =  -0.0051846731;
    double f0_980_img  = -0.26995022;

    Variable v_f0_980_Mass("f0_980_MASS"+Charge, f0_980_MASS);
    Variable v_f0_980_GPP("f0_980_GPP"+Charge, f0_980_GPP);
    Variable v_f0_980_GKK("f0_980_GKK"+Charge, f0_980_GKK);
    Variable v_f0_980_real("f0_980_REAL"+Charge, f0_980_re, 0.01, 0, 0);
    Variable v_f0_980_img("f0_980_IMAG"+Charge, f0_980_img, 0.01, 0, 0);

    auto f0_980 = new Resonances::FLATTE(
        "f0_980"+Charge, v_f0_980_real, v_f0_980_img, v_f0_980_Mass, v_f0_980_GPP, v_f0_980_GKK, 0, PAIR_13 ,  true);
    //P-wave

    double rho770_MASS  = 0.77526;
    double rho770_WIDTH = 0.1478;
    double rho770_amp   = 1.;
    double rho770_img   = 0.;

    Variable v_rho770_Mass("rho770_MASS"+Charge, rho770_MASS);
    Variable v_rho770_Width("rho770_WIDTH"+Charge, rho770_WIDTH);
    Variable v_rho770_real("rho770_REAL"+Charge, rho770_amp);
    Variable v_rho770_img("rho770_IMAG"+Charge, rho770_img);

    auto rho770 = new Resonances::RBW(
        "rho770"+Charge, v_rho770_real, v_rho770_img, v_rho770_Mass, v_rho770_Width, 1, PAIR_13, true, true);

    double rho1450_MASS  = 1.465;
    double rho1450_WIDTH = 0.4;
    double rho1450_amp   = charge_plus ? 0.37*cos(1.99) : 0.37*cos(1.99+0.7);//-0.15060229;
    double rho1450_img   = charge_plus ? 0.37*sin(1.99) : 0.37*sin(1.99+0.7);//0.33796294;

    Variable v_rho1450_Mass("rho1450_MASS"+Charge, rho1450_MASS);
    Variable v_rho1450_Width("rho1450_WIDTH"+Charge, rho1450_WIDTH);
    Variable v_rho1450_real("rho1450_REAL"+Charge, rho1450_amp,0.01,0,0);
    Variable v_rho1450_img("rho1450_IMAG"+Charge, rho1450_img,0.01,0,0);


    auto rho1450 = new Resonances::RBW(
        "rho1450"+Charge, v_rho1450_real, v_rho1450_img, v_rho1450_Mass, v_rho1450_Width, 1, PAIR_13, true, true);

    // D-wave
    double f2_1270_MASS  = 1.2751;
    double f2_1270_WIDTH = 0.1851;
    double f2_1270_amp   = 0.095300878;
    double f2_1270_img   = 0.52136143;

    Variable v_f2_1270_Mass("f2_1270_MASS"+Charge, f2_1270_MASS);
    Variable v_f2_1270_Width("f2_1270_WIDTH"+Charge, f2_1270_WIDTH);
    Variable v_f2_1270_real("f2_1270_REAL"+Charge, f2_1270_amp,0.001,0,0);
    Variable v_f2_1270_img("f2_1270_IMAG"+Charge, f2_1270_img,0.001,0,0);

    auto f2_1270 = new Resonances::RBW(
        "f2"+Charge, v_f2_1270_real, v_f2_1270_img, v_f2_1270_Mass, v_f2_1270_Width, 2, PAIR_13, true, true);

    //swave
    //auto swave = loadPWAResonance("PWAFile.bin");

    
    std::vector<ResonancePdf *> vec_resonances;
    vec_resonances.push_back(rho770);
    vec_resonances.push_back(rho1450);
    vec_resonances.push_back(f0_980);
    vec_resonances.push_back(f2_1270);
    vec_resonances.push_back(NR);

    //vec_resonances.push_back(swave);

    D2KKK.resonances = vec_resonances;

    return new Amp3BodySqDP("signalPDF"+Charge, mprime, thetaprime, eventNumber, D2KKK, eff);
}

void genToyMCB(size_t Nevents, std::string toyName, Amp3BodySqDP *signal, bool charge_pos=true){
    
            auto name= fmt::format("{0}",toyName);
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
            //auto b_jac = t->Branch("jac",&_jac,"jac/D");
            
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
                            //_jac = calc_SqDp_Jacobian(mp, th, Decay_MASS, d1_MASS, d2_MASS, d3_MASS);
                            t->Fill();   
                            accepted++;
                        }
                }
            }

            t->Write("",TObject::kOverwrite);
            f->Close();
            std::cout << "------------------------------------------" << std::endl;
            if(charge_pos)
                std::cout << "toyMC Plus --> " << name.c_str() << " was saved!" << std::endl;
            else
                std::cout << "toyMC Minus --> " << name.c_str() << " was saved!" << std::endl;
            std::cout << "------------------------------------------" << std::endl;

            auto frac = signal->fit_fractions(true);
}

void getData(std::string toyFileName, GooFit::Application &app, DataSet &data, bool charge_plus=true) {
    // load data in a GooFit::dataset

    auto obs = data.getObservables();
    Observable mprime         = obs.at(0);
    Observable thetaprime         = obs.at(1);
    Observable category = obs.at(2);
    Observable eventNumber = obs.at(3);

    auto openRoot = new TFile(toyFileName.c_str());
    auto tree     = (TTree *)openRoot->Get("genResults");
    auto mprime_val(0.);
    auto thetaprime_val(0.);
    auto flags_val(0.);

    tree->SetBranchAddress("mPrime", &mprime_val);
    tree->SetBranchAddress("thPrime", &thetaprime_val);


    for(size_t i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        mprime.setValue(mprime_val);
        thetaprime.setValue(thetaprime_val);
        eventNumber.setValue(data.getNumEvents());
        category.setValue(charge_plus? 0 : 1);
        data.addEvent();
    }
}

void mergeDatasets(UnbinnedDataSet &data1, UnbinnedDataSet &data2, UnbinnedDataSet &data_out){
    // merge two datasets
    auto data1_M = data1.to_matrix<MatrixXd>();
    auto data2_M = data2.to_matrix<MatrixXd>();

    MatrixXd concatenated_matrix(data1_M.rows(), data1_M.cols() + data2_M.cols());
    concatenated_matrix << data1_M, data2_M;

    data_out.from_matrix(concatenated_matrix);

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
    bool randomize              = false;
    size_t Nevents              = 100000;


    auto fit = app.add_subcommand("fit", "fit data");
    fit->add_option("-f,--file", input_data_name, "name_of_file.root");
    fit->add_option("-t,--isToy", is_toy, "Get toyData for fit");
    fit->add_option("-s,--saveToy", save_toy, "save toy in root file");
    fit->add_option("-n,--fitName", fit_name, "name of this fit(useful to save results)");
    fit->add_option("-a,--acc", acc_file, "name of acc file");
    fit->add_option("-b,--bkg", bkg_file, "name of bkg file");
    fit->add_option("-r,--randomize", randomize, "randomize pars");

    auto makeToyMCB = app.add_subcommand("makeToyMCB", "make toy using MC Booster");
    makeToyMCB->add_option("-e,--nevents", Nevents, "number of events");
    makeToyMCB->add_option("-s,--saveToy", save_toy, "save toy in root file");
    makeToyMCB->add_option("-a,--acc", acc_file, "name of acc file");
    makeToyMCB->add_option("-b,--bkg", bkg_file, "name of bkg file");

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

    if(*fit){
        std::cout << "START Fitting data!" << std::endl;

        //Positive PDF
        auto efficiency_plus  = polyEff(mprime, thetaprime);
        auto signal_plus = makesignalpdf(mprime, thetaprime, eventNumber, efficiency_plus);
        auto totalpdf_plus    = new ProdPdf("totalpdf_Plus", {signal_plus});

        //Negative PDF
        auto efficiency_minus  = polyEff(mprime, thetaprime, false);
        auto signal_minus = makesignalpdf(mprime, thetaprime, eventNumber, efficiency_minus, false);
        auto totalpdf_minus   = new ProdPdf("totalpdf_Minus", {signal_minus});

        auto stepFunction =  BinTransformPdf("stepFunction", {category}, {-0.5}, {1}, {2});
       
        std::vector<GooPdf *> comps;
        comps.push_back(totalpdf_plus);
        comps.push_back(totalpdf_minus);
        auto finalPDF_comb = MappedPdf("finalPDF_comb", &stepFunction, comps);

        //load datasets
        UnbinnedDataSet dataPlus({mprime, thetaprime, category, eventNumber});
        UnbinnedDataSet dataMinus({mprime, thetaprime, category, eventNumber});
        UnbinnedDataSet dataMerged({mprime, thetaprime, category, eventNumber});
        getData("MC/MC_Plus.root", app, dataPlus);
        printf("Number of events in dataset plus: %d \n",dataPlus.getNumEvents());
        getData("MC/MC_Minus.root", app, dataMinus, false);
        printf("Number of events in dataset minus: %d \n",dataMinus.getNumEvents());
        mergeDatasets(dataPlus, dataMinus, dataMerged); 
        printf("Number of events in combined dataset: %d \n",dataMerged.getNumEvents());

        signal_plus->setDataSize(dataPlus.getNumEvents(),4,0);
        signal_minus->setDataSize(dataMinus.getNumEvents(),4,dataPlus.getNumEvents());

        finalPDF_comb.setData(&dataMerged);

        if(randomize){
            std::random_device rd; 
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> dis(0.,1.);
            auto parameters = finalPDF_comb.getParameters();
            for(auto p: parameters){
                if(p.IsFixed())
                    continue;
                p.setValue(dis(gen));
            }
        }

        FitManager datapdf(&finalPDF_comb);
        datapdf.setVerbosity(2);
        datapdf.setMaxCalls(200000);
        datapdf.setStrategy(2);
        // Start fit
        auto func_min = datapdf.fit();


        {
            signal_plus->resetCacheCounter();
            std::cout << "Fit fractions plus: " << std::endl;
            //signal_plus->fit_fractions(true);
            toyName = fmt::format("Fit/{0}/MC_Plus.root", fit_name);
            genToyMCB(Nevents, toyName, signal_plus, true);
        }
        std::cout << "------------------------------------------" << std::endl;
        {
            std::cout << "Fit fractions minus: " << std::endl;
            //signal_minus->fit_fractions(true);
            toyName = fmt::format("Fit/{0}/MC_Minus.root", fit_name);
            genToyMCB(Nevents, toyName, signal_minus, false);
        }

    }

    if(*makeToyMCB) {
        
        {
            //Positive PDF
            auto efficiency_plus  = polyEff(mprime, thetaprime);
            auto signal_plus = makesignalpdf(mprime, thetaprime, eventNumber, efficiency_plus);
            auto totalpdf_plus    = new ProdPdf("totalpdf_Plus", {signal_plus});
           
            toyName = "MC/MC_Plus.root";
            genToyMCB(Nevents, toyName, signal_plus, true);
            signal_plus->resetCacheCounter();
        }

        {
            //Negative PDF
            auto efficiency_minus  = polyEff(mprime, thetaprime, false);
            auto signal_minus = makesignalpdf(mprime, thetaprime, eventNumber, efficiency_minus, false);
            auto totalpdf_minus   = new ProdPdf("totalpdf_Minus", {signal_minus});

            toyName = "MC/MC_Minus.root";
            genToyMCB(Nevents, toyName, signal_minus, false);
            signal_minus->resetCacheCounter();
        }

    }




}