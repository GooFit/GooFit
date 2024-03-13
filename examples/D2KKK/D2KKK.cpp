// Dalitz plot analysis of the decay D± → K±K±K∓
// arXiv:1902.05884v3
//  ROOT stuff
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
#include <goofit/PDFs/physics/Amp3Body.h>
#include <goofit/PDFs/physics/DalitzVetoPdf.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>
#include <goofit/PDFs/physics/DalitzPlotter.h>
#include <goofit/FunctorWriter.h>

// Matrix
#include <Eigen/Core>
#include <Eigen/LU>

#define torad(x) (x * M_PI / 180.)

using namespace std;
using namespace GooFit;
using namespace ROOT;

// D_MASS (PDG)
double D_MASS = 1.86966;
double k_MASS = 0.493677;

double Decay_MASS = D_MASS;
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
Observable s13("s13", s12_min, s12_max);
Observable s23("s23", s13_min, s13_max);
Observable s12("s12", s12_min, s12_max);
EventNumber eventNumber("eventNumber");

GooPdf *polyEff(Observable s13, Observable s23) {
    vector<Variable> offsets;
    vector<Observable> observables;
    vector<Variable> coefficients;
    Variable constantOne("c1", 1);
    Variable constantZero("c0", 0);
    observables.push_back(s13);
    observables.push_back(s23);
    offsets.push_back(constantZero);
    offsets.push_back(constantZero);
    coefficients.push_back(constantOne);
    PolynomialPdf *eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0); // No efficiency

    return eff;
}

Amp3Body *makesignalpdf(Observable s13, Observable s23, EventNumber eventNumber, GooPdf *eff = 0, bool doRandom=false) {
    // set up the decay channel
    DecayInfo3 D2KKK;
    D2KKK.motherMass          = D_MASS;
    D2KKK.daug1Mass           = d1_MASS;
    D2KKK.daug2Mass           = d2_MASS;
    D2KKK.daug3Mass           = d3_MASS;
    D2KKK.meson_radius        = 4.5; // GeV^-1
    D2KKK.mother_meson_radius = 4.5; // GeV^-1

    std::random_device rd; 
    std::mt19937 gen(rd());
    auto randValue = std::uniform_real_distribution<double>(-2.,+2.);
    auto randMassWidth = std::uniform_real_distribution<double>(0.,1.);

    // parameters from Laura++
    double f0_980_MASS = 0.965;
    double f0_980_GPP  = 0.165;
    double f0_980_GKK  = 4.21;
    double f0_980_re   = sqrt(6.17 / 23.7) * 3.12 * cos(-58.9 * M_PI / 180.);
    double f0_980_img  = sqrt(6.17 / 23.7) * 3.12 * sin(-58.9 * M_PI / 180.);

    std::cout << "f0_980_re = " << f0_980_re << std::endl;
    std::cout << "f0_980_img = " << f0_980_img << std::endl;
    std::cout << "f0_980_MASS = " << f0_980_MASS << std::endl;
    std::cout << "f0_980_GPP = " << f0_980_GPP << std::endl;
    std::cout << "f0_980_GKK = " << f0_980_GKK << std::endl;

    if(doRandom){
        f0_980_re = f0_980_re + randValue(gen);
        f0_980_img = f0_980_img + randValue(gen);
    }


   

    Variable v_f0_980_Mass("f0_980_MASS", f0_980_MASS);
    Variable v_f0_980_GPP("f0_980_GPP", f0_980_GPP);
    Variable v_f0_980_GKK("f0_980_GKK", f0_980_GKK);
    Variable v_f0_980_real("f0_980_REAL", f0_980_re, 0.01, 0, 0);
    Variable v_f0_980_img("f0_980_IMAG", f0_980_img, 0.01, 0, 0);

    auto f0_980 = new Resonances::FLATTE(
        "f0_980", v_f0_980_real, v_f0_980_img, v_f0_980_Mass, v_f0_980_GPP, v_f0_980_GKK, 0, PAIR_13 ,  true);

    double f0_1370_MASS  = 1.422;
    double f0_1370_WIDTH = 0.324;
    double f0_1370_amp   = sqrt(6.17 / 25.4) * 3.46 * cos(13.1 * M_PI / 180.);
    double f0_1370_img   = sqrt(6.17 / 25.4) * 3.46 * sin(13.1 * M_PI / 180.);
    std::cout << "f0_1370_re = " << f0_1370_amp << std::endl;
    std::cout << "f0_1370_img = " << f0_1370_img << std::endl;
    std::cout << "f0_1370_MASS = " << f0_1370_MASS << std::endl;
    std::cout << "f0_1370_WIDTH = " << f0_1370_WIDTH << std::endl;

    if(doRandom){
        f0_1370_amp = f0_1370_amp + randValue(gen);
        f0_1370_img = f0_1370_img + randValue(gen);
        f0_1370_MASS = f0_1370_MASS*randMassWidth(gen);
        f0_1370_WIDTH = f0_1370_WIDTH*randMassWidth(gen);
    }

    Variable v_f0_1370_Mass("f0_1370_MASS", f0_1370_MASS, 0.01, 0, 0);
    Variable v_f0_1370_Width("f0_1370_WIDTH", f0_1370_WIDTH, 0.01, 0, 0);
    Variable v_f0_1370_real("f0_1370_REAL", f0_1370_amp, 0.01, 0, 0);
    Variable v_f0_1370_img("f0_1370_IMAG", f0_1370_img, 0.01, 0, 0);

    auto f0_1370 = new Resonances::RBW(
        "f0_1370", v_f0_1370_real, v_f0_1370_img, v_f0_1370_Mass, v_f0_1370_Width, 0, PAIR_13, true, true);

    double phi1020_MASS  = 1.019461;
    double phi1020_WIDTH = 0.00429;
    double phi1020_amp   = 1.;
    double phi1020_img   = 0.;

    Variable v_phi1020_Mass("phi1020_MASS", phi1020_MASS);
    Variable v_phi1020_Width("phi1020_WIDTH", phi1020_WIDTH);
    Variable v_phi1020_real("phi1020_REAL", phi1020_amp);
    Variable v_phi1020_img("phi1020_IMAG", phi1020_img);

    auto phi1020 = new Resonances::RBW(
        "phi1020", v_phi1020_real, v_phi1020_img, v_phi1020_Mass, v_phi1020_Width, 1, PAIR_13, true, true);

    // If you want include a resonance in your model, just push into the vector 'vec_resonances'
    std::vector<ResonancePdf *> vec_resonances;

    vec_resonances.push_back(f0_980);
    vec_resonances.push_back(f0_1370);
    vec_resonances.push_back(phi1020);

    D2KKK.resonances = vec_resonances;

    return new Amp3Body("signalPDF", s13, s23, eventNumber, D2KKK, eff);
}

void getData(std::string toyFileName, GooFit::Application &app, DataSet &data, bool toy) {
    // load data in a GooFit::dataset

    toyFileName = app.get_filename(toyFileName, "MC/");

    auto obs               = data.getObservables();
    Observable s13         = obs.at(0);
    Observable s23         = obs.at(1);
    Observable eventNumber = obs.at(2);

    auto openRoot = new TFile(toyFileName.c_str());
    auto tree     = (TTree *)openRoot->Get("DecayTree");
    auto s13_val(0.);
    auto s23_val(0.);

    if(toy) {
        tree->SetBranchAddress("s13", &s13_val);
        tree->SetBranchAddress("s23", &s23_val);
    } else {
        tree->SetBranchAddress("s13_pipi_DTF", &s13_val);
        tree->SetBranchAddress("s23_pipi_DTF", &s23_val);
    }

    size_t j = 0;
    for(size_t i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        s13.setValue(s13_val);
        s23.setValue(s23_val);
        eventNumber.setValue(data.getNumEvents());
        if((s13.getValue() < s13.getUpperLimit()) && (s23.getValue() < s23.getUpperLimit())
           && (s13.getValue() > s13.getLowerLimit()) && (s23.getValue() > s23.getLowerLimit())) {
            data.addEvent();
            if(j < 10)
                printf("[%d] = (%f , %f)\n", i, s13.getValue(), s23.getValue());
            j++;
        }
    }
}

void to_root(UnbinnedDataSet &toyMC, std::string name) {
    // save GooFit::Dataset in a root file
    auto obs               = toyMC.getObservables();
    Observable s13         = obs.at(0);
    Observable s23         = obs.at(1);
    Observable eventNumber = obs.at(2);

    double _s12, _s13, _s23;
    auto f     = new TFile(name.c_str(), "recreate");
    auto t     = new TTree("DecayTree", "toyMC");
    auto b_s12 = t->Branch("s12", &_s12, "s12/D");
    auto b_s13 = t->Branch("s13", &_s13, "s13/D");
    auto b_s23 = t->Branch("s23", &_s23, "s23/D");

    for(int i = 0; i < toyMC.getNumEvents(); i++) {
        toyMC.loadEvent(i);
        t->GetEntry(i);
        _s13 = s13.getValue();
        _s23 = s23.getValue();
        _s12 = POW2(Decay_MASS) + POW2(d1_MASS) + POW2(d2_MASS) + POW2(d3_MASS) - s13.getValue() - s23.getValue();
        t->Fill();
    }
    t->Write("", TObject::kOverwrite);
    f->Write("", TObject::kOverwrite);
    f->Close();
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "toyMC --> " << name.c_str() << " was saved!" << std::endl;
    std::cout << "------------------------------------------" << std::endl;
}

Amp3Body *runFit(GooPdf *totalPdf, Amp3Body *signal, UnbinnedDataSet *data, std::string name) {
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
    size_t Nevents              = 1000000;
    bool randpars = false;

    auto fit = app.add_subcommand("fit", "fit data");
    fit->add_option("-f,--file", input_data_name, "name_of_file.root");
    fit->add_option("-t,--isToy", is_toy, "Get toyData for fit");
    fit->add_option("-s,--saveToy", save_toy, "save toy in root file");
    fit->add_option("-n,--fitName", fit_name, "name of this fit(useful to save results)")->required(true);
    fit->add_option("-a,--acc", acc_file, "name of acc file");
    fit->add_option("-b,--bkg", bkg_file, "name of bkg file");
    fit->add_option("-d,--disable-acc-bkg", no_acc_and_bkg, "disable-acc-bkg");
    fit->add_option("-r,--randomize", randpars, "randomize parameters");

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

    s13.setNumBins(bins);
    s23.setNumBins(bins);

    auto efficiency  = polyEff(s13, s23);
    auto background  = nullptr;
    GooPdf *totalpdf = nullptr;

    // to be implemented...
    /*if(!no_acc_and_bkg){
         efficiency = makeHistogramPdf(efffile,effhist,s13,s23,true,false,false);
         background = makeHistogramPdf(bkgfile,bkghist,s13,s23,false,false,false);
         auto signal = makesignalpdf(s13, s23, eventNumber,efficiency);
         totalpdf = new AddPdf("totalpdf", Variable("frac",0.93), signal, background) ;

    }*/

    auto signal = makesignalpdf(s13, s23, eventNumber, efficiency, randpars);
    totalpdf    = new ProdPdf("totalpdf", {signal});

    if(*makeToy) {
        DalitzPlotter dplotter{totalpdf, signal};
        UnbinnedDataSet data({s13, s23, eventNumber});
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

        std::cout << "Fit Fractions Interference" << '\n';

        // signal->normalize();
        auto frac = signal->fit_fractions(true);

        return 0;
    }

    /* it uses accept and reject method (very inefficient for B phase space)
    if(*makeToy) {
        auto tuple = signal->GenerateSig(Nevents,200323);
        if(save_toy) {
            auto name= fmt::format("MC/{0}",toyName);
            double _s12, _s13,_s23;

            auto f = new TFile(name.c_str(),"recreate");
            auto t = new TTree("DecayTree","toyMC");
            auto b_s12 = t->Branch("s13",&_s12,"s13/D");
            auto b_s13 = t->Branch("s23",&_s13,"s23/D");
            auto b_s23 = t->Branch("s23",&_s23,"s23/D");

            auto variables = std::get<1>(tuple);
            auto weights   = std::get<2>(tuple);
            auto flags     = std::get<3>(tuple);

            for(int i = 0; i < weights.size(); i++){
                if(flags[i] == 1) {
                    t->GetEntry(i);
                    _s12 = (*(variables[0]))[i];
                    _s13 = (*(variables[2]))[i];
                    _s23 = (*(variables[1]))[i];
                    t->Fill();
                }
            }
            t->Write("",TObject::kOverwrite);
            f->Write();
            f->Close();
            std::cout << "------------------------------------------" << std::endl;
            std::cout << "toyMC --> " << name.c_str() << " was saved!" << std::endl;
            std::cout << "------------------------------------------" << std::endl;
        }
        return 0;
    }*/

    if(*fit) {
        auto command = fmt::format("mkdir -p Fit/{0}", fit_name);
        if(system(command.c_str()) != 0)
            throw GooFit::GeneralError("Making directory failed");

        std::cout << "------------------------------------------" << std::endl;
        std::cout << "Reading file --> " << input_data_name << std::endl;
        std::cout << "------------------------------------------" << std::endl;
        UnbinnedDataSet data({s13, s23, eventNumber});
        getData(input_data_name, app, data, is_toy);
        std::cout << "------------------------------------------" << std::endl;
        std::cout << "Num Entries Loaded =  " << data.getNumEvents() << std::endl;
        std::cout << "------------------------------------------" << std::endl;
        auto output_signal = runFit(totalpdf, signal, &data, fit_name);

        std::cout << "norm = " << output_signal->normalize() << std::endl;

        auto plotpath = fmt::format("Fit/{0}/plots.root", fit_name);
        DalitzPlotter dplotter{totalpdf, signal};
        dplotter.Plot(plotpath.c_str(), &data, 120);
    }
}
