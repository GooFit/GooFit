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
double d1_MASS    = k_MASS;
double d2_MASS    = k_MASS;
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

Amp3BodySqDP *makesignalpdf(Observable mprime, Observable thetaprime, EventNumber eventNumber, GooPdf *eff = 0) {
    // set up the decay channel
    DecayInfo3 D2KKK;
    D2KKK.motherMass          = B_MASS;
    D2KKK.daug1Mass           = d1_MASS;
    D2KKK.daug2Mass           = d2_MASS;
    D2KKK.daug3Mass           = d3_MASS;
    D2KKK.meson_radius        = 4.0; // GeV^-1
    D2KKK.mother_meson_radius = 5.0; // GeV^-1


    double f2p_1525_MASS  = 1.5174;
    double f2p_1525_WIDTH = 0.086;
    double f2p_1525_amp   = 0.;
    double f2p_1525_img   = 1.;

    Variable v_f2p_1525_Mass("f2p_1525_MASS", f2p_1525_MASS);
    Variable v_f2p_1525_Width("f2p_1525_WIDTH", f2p_1525_WIDTH);
    Variable v_f2p_1525_real("f2p_1525_REAL", f2p_1525_amp, 0.01, 0, 0);
    Variable v_f2p_1525_img("f2p_1525_IMAG", f2p_1525_img, 0.01, 0, 0);

    auto f2p_1525 = new Resonances::RBW(
        "f2p_1525", v_f2p_1525_real, v_f2p_1525_img, v_f2p_1525_Mass, v_f2p_1525_Width, 2, PAIR_13, true, false);

    double phi1020_MASS  = 1.019461;
    double phi1020_WIDTH = 0.00429;
    double phi1020_amp   = 1.;
    double phi1020_img   = 0.;

    Variable v_phi1020_Mass("phi1020_MASS", phi1020_MASS);
    Variable v_phi1020_Width("phi1020_WIDTH", phi1020_WIDTH);
    Variable v_phi1020_real("phi1020_REAL", phi1020_amp);
    Variable v_phi1020_img("phi1020_IMAG", phi1020_img);

    auto phi1020 = new Resonances::RBW(
        "phi1020", v_phi1020_real, v_phi1020_img, v_phi1020_Mass, v_phi1020_Width, 1, PAIR_13, true, false);

    // If you want include a resonance in your model, just push into the vector 'vec_resonances'

    auto nonres = new Resonances::NonRes("NonRes",Variable("re",1.,0.01,0,0),Variable("im",0.,0.01,0,0));

    std::vector<ResonancePdf *> vec_resonances;

    vec_resonances.push_back(phi1020);
    vec_resonances.push_back(f2p_1525);
    vec_resonances.push_back(nonres);

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

    printf("NEntries = %f \n",tree->GetEntries());

    tree->SetBranchAddress("mPrime", &mprime_val);
    tree->SetBranchAddress("thPrime", &thetaprime_val);
   
    size_t j = 0;
    for(size_t i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
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
        if(_thetaprime>0.5)
            _thetaprime = 1.0 -_thetaprime;
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
    size_t Nevents              = 100000;

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
    GooPdf *totalpdf = nullptr;

    auto signal = makesignalpdf(mprime, thetaprime, eventNumber, efficiency);
    //signal->setNumNormEvents(1000000);
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