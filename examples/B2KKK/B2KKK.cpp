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
#include <goofit/PDFs/basic/StepPdf.h>
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

using namespace std;
using namespace GooFit;
using namespace ROOT;

// B_MASS (PDG)
double B_MASS = 5.27934;
double k_MASS = 0.493677;

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
Observable s12("s12", s12_min, s12_max);
Observable s13("s13", s13_min, s13_max);
Observable s23("s23", s23_min, s23_max);
EventNumber eventNumber("eventNumber");

GooPdf *polyEff() {
    vector<Variable> offsets;
    vector<Observable> observables;
    vector<Variable> coefficients;
    Variable constantOne("c1", 1);
    Variable constantZero("c0", 0);
    observables.push_back(s12);
    observables.push_back(s13);
    offsets.push_back(constantZero);
    offsets.push_back(constantZero);
    coefficients.push_back(constantOne);
    PolynomialPdf *eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0); // No efficiency

    return eff;
}

GooPdf *vetoPdf() {
    // for Dzero veto
    const fptype D0Mass = 1.86484;
    const fptype range  = 0.030;
    fptype veto_min     = D0Mass * D0Mass;
    fptype veto_max     = s12_max;
    VetoInfo D0VetoInfoLow(Variable("veto_Dzero_low_min", veto_min), Variable("veto_Dzero_low_max", veto_max), PAIR_12);
    VetoInfo D0VetoInfoHigh(
        Variable("veto_Dzero_high_min", veto_min), Variable("veto_Dzero_high_max", veto_max), PAIR_13);

    vector<VetoInfo> vetos;
    vetos.push_back(D0VetoInfoLow);
    vetos.push_back(D0VetoInfoHigh);

    auto ret
        = new DalitzVetoPdf("vetoPdf", s12, s13, Mother_Mass, Daughter1_Mass, Daughter2_Mass, Daughter3_Mass, vetos);

    return ret;
}

// https://arxiv.org/pdf/hep-ex/0412066.pdf
Amp3Body *makesignalpdf(GooPdf *eff = 0) {
    // set up the decay channel
    DecayInfo3 B2KKK;
    B2KKK.motherMass          = B_MASS;
    B2KKK.daug1Mass           = d1_MASS;
    B2KKK.daug2Mass           = d2_MASS;
    B2KKK.daug3Mass           = d3_MASS;
    B2KKK.meson_radius        = 1.5; // GeV^-1
    B2KKK.mother_meson_radius = 5.0;

    double phi1020_MASS  = 1.019461;
    double phi1020_WIDTH = 0.00429;
    double phi1020_real  = 0.15 * cos(-200 * M_PI / 180.);
    double phi1020_img   = 0.15 * sin(-200 * M_PI / 180.);

    Variable v_phi1020_Mass("phi1020_MASS", phi1020_MASS);
    Variable v_phi1020_Width("phi1020_WIDTH", phi1020_WIDTH);
    Variable v_phi1020_real("phi1020_REAL", phi1020_real, 0.01, 0, 0);
    Variable v_phi1020_img("phi1020_IMAG", phi1020_img, 0.01, 0, 0);

    auto phi1020 = new Resonances::RBW(
        "phi1020", v_phi1020_real, v_phi1020_img, v_phi1020_Mass, v_phi1020_Width, 1, PAIR_12, true, false);

    double rho1450_MASS  = 1.465;
    double rho1450_WIDTH = 0.4;
    double rho1450_real  = -0.573;
    double rho1450_img   = -1.139;

    Variable v_rho1450_Mass("rho1450_MASS", rho1450_MASS);
    Variable v_rho1450_Width("rho1450_WIDTH", rho1450_WIDTH);
    Variable v_rho1450_real("rho1450_REAL", rho1450_real, 0.01, 0, 0);
    Variable v_rho1450_img("rho1450_IMAG", rho1450_img, 0.01, 0, 0);

    auto rho1450 = new Resonances::RBW(
        "rho1450", v_rho1450_real, v_rho1450_img, v_rho1450_Mass, v_rho1450_Width, 1, PAIR_12, true, false);

    double rho1700_MASS  = 1.720;
    double rho1700_WIDTH = 0.250;
    double rho1700_real  = 0.449;
    double rho1700_img   = -0.743;

    Variable v_rho1700_Mass("rho1700_MASS", rho1700_MASS);
    Variable v_rho1700_Width("rho1700_WIDTH", rho1700_WIDTH);
    Variable v_rho1700_real("rho1700_REAL", rho1700_real, 0.01, 0, 0);
    Variable v_rho1700_img("rho1700_IMAG", rho1700_img, 0.01, 0, 0);

    auto rho1700 = new Resonances::RBW(
        "rho1700", v_rho1700_real, v_rho1700_img, v_rho1700_Mass, v_rho1700_Width, 1, PAIR_12, true, false);

    double f0980_MASS = 0.965;
    double f0980_gpi  = 0.165;
    double f0980_gk   = 4.21 * f0980_gpi;
    double f0980_real = 1.;
    double f0980_img  = 0.;

    Variable v_f0980_Mass("f0980_MASS", f0980_MASS);
    Variable v_f0980_gpi("f0980_WIDTH", f0980_gpi);
    Variable v_f0980_gk("f0980_WIDTH", f0980_gk);
    Variable v_f0980_real("f0980_REAL", f0980_real, 0.01, 0, 0);
    Variable v_f0980_img("f0980_IMAG", f0980_img, 0.01, 0, 0);

    auto f0980 = new Resonances::FLATTE(
        "f0980", v_f0980_real, v_f0980_img, v_f0980_Mass, v_f0980_gpi, v_f0980_gk, PAIR_12, true);

    double f01500_MASS  = 1.491;
    double f01500_WIDTH = 0.145;
    double f01500_real  = 1.;
    double f01500_img   = 0.;

    Variable v_f01500_Mass("f01500_MASS", f01500_MASS);
    Variable v_f01500_Width("f01500_WIDTH", f01500_WIDTH);
    Variable v_f01500_real("f01500_REAL", f01500_real, 0.01, 0, 0);
    Variable v_f01500_img("f01500_IMAG", f01500_img, 0.01, 0, 0);

    auto f01500 = new Resonances::RBW(
        "f01500", v_f01500_real, v_f01500_img, v_f01500_Mass, v_f01500_Width, 0, PAIR_12, true, false);

    fptype Xc0_I_MASS  = 3.41475;
    fptype Xc0_I_WIDTH = 0.0108;
    fptype Xc0_I_amp   = 0.08 * cos(127. * M_PI / 180.);
    fptype Xc0_I_img   = 0.08 * sin(127. * M_PI / 180.);
    Variable v_Xc0_I_Mass("Xc0_I_MASS", Xc0_I_MASS);
    Variable v_Xc0_I_Width("Xc0_I_WIDTH", Xc0_I_WIDTH);
    Variable v_Xc0_I_real("Xc0_I_REAL", Xc0_I_amp);
    Variable v_Xc0_I_img("Xc0_I_IMAG", Xc0_I_img);

    auto Xc0_I
        = new Resonances::RBW("Xc0_I", v_Xc0_I_real, v_Xc0_I_img, v_Xc0_I_Mass, v_Xc0_I_Width, 0, PAIR_12, true, false);

    Variable v_NR_real("NR_REAL", 0.65 * cos(-61 * M_PI / 180), 0.01, 0, 0);
    Variable v_NR_img("NR_IMAG", 0.65 * sin(-61 * M_PI / 180), 0.01, 0, 0);
    Variable v_NR_alpha("NR_alpha", 0.116);

    auto NR = new Resonances::NonRes3k("NR", v_NR_real, v_NR_img, v_NR_alpha);

    // If you want include a resonance in your model, just push into the vector 'vec_resonances'
    std::vector<ResonancePdf *> vec_resonances;

    vec_resonances.push_back(phi1020);
    vec_resonances.push_back(rho1450);
    vec_resonances.push_back(rho1700);
    vec_resonances.push_back(f0980);
    vec_resonances.push_back(f01500);
    vec_resonances.push_back(Xc0_I);
    vec_resonances.push_back(NR);

    B2KKK.resonances = vec_resonances;

    return new Amp3Body("signalPDF", s12, s13, eventNumber, B2KKK, eff);
}

void getData(DataSet &data) {
    // load data in a GooFit::dataset

    auto openRoot = new TFile("B2KKK.root");
    auto tree     = (TTree *)openRoot->Get("DecayTree");
    auto s12_val(0.);
    auto s13_val(0.);

    tree->SetBranchAddress("s12", &s12_val);
    tree->SetBranchAddress("s13", &s13_val);

    size_t j = 0;
    for(int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        s12.setValue(s12_val);
        s13.setValue(s13_val);
        eventNumber.setValue(data.getNumEvents());
        if((s12.getValue() < s12.getUpperLimit()) && (s13.getValue() < s13.getUpperLimit())
           && (s12.getValue() > s12.getLowerLimit()) && (s13.getValue() > s13.getLowerLimit())) {
            data.addEvent();
            if(j < 10)
                printf("[%d] = (%f , %f)\n", i, s12.getValue(), s13.getValue());
            j++;
        }
    }
}

void to_root(UnbinnedDataSet &toyMC, std::string name) {
    // save GooFit::Dataset in a root file
    auto obs               = toyMC.getObservables();
    Observable s12         = obs.at(0);
    Observable s13         = obs.at(1);
    Observable eventNumber = obs.at(2);

    double _s12, _s13, _s23;
    auto f = new TFile(name.c_str(), "recreate");
    auto t = new TTree("DecayTree", "toyMC");
    t->Branch("s12", &_s12, "s12/D");
    t->Branch("s13", &_s13, "s13/D");
    t->Branch("s23", &_s23, "s23/D");

    for(int i = 0; i < toyMC.getNumEvents(); i++) {
        toyMC.loadEvent(i);
        t->GetEntry(i);
        _s12 = s12.getValue();
        _s13 = s13.getValue();
        _s23 = POW2(Decay_MASS) + POW2(d1_MASS) + POW2(d2_MASS) + POW2(d3_MASS) - s12.getValue() - s13.getValue();
        t->Fill();
    }
    t->Write(0, TObject::kOverwrite);
    f->Write();
    f->Close();
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "toyMC --> " << name.c_str() << " was saved!" << std::endl;
    std::cout << "------------------------------------------" << std::endl;
}

void runfit1() {
    auto efficiency = polyEff();
    efficiency->setParameterConstantness(true);
    auto signal   = makesignalpdf(efficiency);
    auto totalpdf = new ProdPdf("totalpdf", {signal});

    UnbinnedDataSet data({s12, s13, eventNumber});
    getData(data);

    totalpdf->setData(&data);
    signal->setDataSize(data.getNumEvents());
    auto fitter = FitManager(totalpdf);
    fitter.setVerbosity(1);
    fitter.setMaxCalls(200000);
    auto func_min = fitter.fit();
    writeToFile(totalpdf, "B2KKK_fitwithoutveto.txt");
}

void runfit2() {
    auto efficiency = polyEff();
    efficiency->setParameterConstantness(true);
    auto vetopdf = vetoPdf();
    vetopdf->setParameterConstantness(true);
    auto effwithveto      = new ProdPdf("effwithveto", {vetopdf, efficiency});
    auto signalwithveto   = makesignalpdf(effwithveto);
    auto totalpdfwithveto = new ProdPdf("totalpdfwithveto", {signalwithveto});

    UnbinnedDataSet data({s12, s13, eventNumber});
    getData(data);

    totalpdfwithveto->setData(&data);
    signalwithveto->setDataSize(data.getNumEvents());
    auto fitter = FitManager(totalpdfwithveto);
    fitter.setVerbosity(1);
    fitter.setMaxCalls(200000);
    auto func_min = fitter.fit();
    writeToFile(totalpdfwithveto, "B2KKK_fitwithveto.txt");
}

UnbinnedDataSet fillData() {
    auto efficiency = polyEff();
    efficiency->setParameterConstantness(true);
    auto signal   = makesignalpdf(efficiency);
    auto totalpdf = new ProdPdf("totalpdf", {signal});

    writeToFile(totalpdf, "B2KKK_initialpars.txt");

    DalitzPlotter dplotter{totalpdf, signal};
    UnbinnedDataSet data({s12, s13, eventNumber});
    dplotter.fillDataSetMC(data, 2000000);

    auto toyhist = (TH2F *)dplotter.make2D();
    toyhist->Rebin2D(2);
    TCanvas foo("foo", "", 1000, 500);
    gStyle->SetOptStat(0);
    foo.Divide(2);
    foo.cd(1);
    gPad->SetLogz();
    toyhist->Draw("colz");
    foo.cd(2);
    gPad->SetLogy();
    auto toys12 = (TH1F *)toyhist->ProjectionX("proj");
    toys12->Draw("hist");
    foo.SaveAs("B2KKK_Simple_toyData.png");

    to_root(data, "B2KKK.root");

    return data;
}

UnbinnedDataSet fillData2() {
    auto efficiency = polyEff();
    efficiency->setParameterConstantness(true);
    auto vetopdf = vetoPdf();
    vetopdf->setParameterConstantness(true);

    auto effwithveto = new ProdPdf("effwithveto", {efficiency, vetopdf});
    auto signal      = makesignalpdf(effwithveto);

    auto totalpdf = new ProdPdf("totalpdf", {signal});

    DalitzPlotter dplotter{totalpdf, signal};
    UnbinnedDataSet data({s12, s13, eventNumber});
    dplotter.fillDataSetMC(data, 2000000);

    auto toyhist = (TH2F *)dplotter.make2D();
    toyhist->Rebin2D(2);
    TCanvas foo("foo", "", 1000, 500);
    gStyle->SetOptStat(0);
    foo.Divide(2);
    foo.cd(1);
    gPad->SetLogz();
    toyhist->Draw("colz");
    foo.cd(2);
    gPad->SetLogy();
    auto toys12 = (TH1F *)toyhist->ProjectionX("proj");
    toys12->Draw("hist");
    foo.SaveAs("B2KKK_Simple_toyData_withveto.png");

    to_root(data, "B2KKK_Veto.root");

    return data;
}

int main(int argc, char **argv) {
    GooFit::Application app{"B2KKK", argc, argv};
    auto makeToy     = app.add_subcommand("makeToy", "make a toy");
    auto makeToyveto = app.add_subcommand("makeToyveto", "make toy wth veto");
    auto _runfit1    = app.add_subcommand("runfit1", "run fit 1");
    auto _runfit2    = app.add_subcommand("runfit2", "run fit 2");
    GOOFIT_PARSE(app);

    s12.setNumBins(bins);
    s13.setNumBins(bins);

    if(*makeToy) {
        fillData();
    }

    if(*makeToyveto) {
        fillData2();
    }

    if(*_runfit1) {
        runfit1();
    }

    if(*_runfit2) {
        runfit2();
    }
}
