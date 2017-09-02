// ROOT stuff
#include "TRandom.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TText.h"
#include "TLine.h"

// System stuff
#include <fstream>
#include <sys/time.h>
#include <sys/times.h>

// GooFit stuff
#include "goofit/Application.h"
#include "goofit/Variable.h"
#include "goofit/PDFs/basic/PolynomialPdf.h"
#include "goofit/PDFs/physics/DalitzPlotPdf.h"
#include "goofit/PDFs/physics/DalitzVetoPdf.h"
#include "goofit/PDFs/physics/ResonancePdf.h"
#include "goofit/PDFs/combine/AddPdf.h"
#include "goofit/PDFs/combine/ProdPdf.h"
#include "goofit/PDFs/GooPdf.h"
#include "goofit/FitManager.h"
#include "goofit/UnbinnedDataSet.h"

using namespace std;
using namespace GooFit;

TCanvas *foo;
TCanvas *foodal;
UnbinnedDataSet *data = 0;

Variable *m12                 = 0;
Variable *m13                 = 0;
CountingVariable *eventNumber = 0;
bool fitMasses                = false;
Variable *fixedRhoMass        = new Variable("rho_mass", 0.7758, 0.01, 0.7, 0.8);
Variable *fixedRhoWidth       = new Variable("rho_width", 0.1503, 0.01, 0.1, 0.2);

const fptype _mD0       = 1.86484;
const fptype _mD02      = _mD0 * _mD0;
const fptype _mD02inv   = 1. / _mD02;
const fptype piPlusMass = 0.13957018;
const fptype piZeroMass = 0.1349766;

// Constants used in more than one PDF component.
Variable *motherM = new Variable("motherM", _mD0);
Variable *chargeM = new Variable("chargeM", piPlusMass);
Variable *neutrlM = new Variable("neutrlM", piZeroMass);
Variable *massSum
    = new Variable("massSum", _mD0 *_mD0 + 2 * piPlusMass * piPlusMass + piZeroMass * piZeroMass); // = 3.53481
Variable *constantOne  = new Variable("constantOne", 1);
Variable *constantZero = new Variable("constantZero", 0);

GooPdf *kzero_veto = 0;

fptype cpuGetM23(fptype massPZ, fptype massPM) {
    return (_mD02 + piZeroMass * piZeroMass + piPlusMass * piPlusMass + piPlusMass * piPlusMass - massPZ - massPM);
}

void getToyData(std::string toyFileName, GooFit::Application &app) {
    toyFileName = app.get_filename(toyFileName, "examples/dalitz");

    TH2F dalitzplot("dalitzplot",
                    "",
                    m12->getNumBins(),
                    m12->getLowerLimit(),
                    m12->getUpperLimit(),
                    m13->getNumBins(),
                    m13->getLowerLimit(),
                    m13->getUpperLimit());
    std::vector<Variable *> vars;
    vars.push_back(m12);
    vars.push_back(m13);
    vars.push_back(eventNumber);
    data = new UnbinnedDataSet(vars);

    std::ifstream reader(toyFileName);
    std::string buffer;

    while(!reader.eof()) {
        reader >> buffer;

        if(buffer == "====")
            break;

        std::cout << buffer;
    }

    double dummy = 0;

    while(!reader.eof()) {
        reader >> dummy;
        reader >> dummy; // m23, m(pi+ pi-), called m12 in processToyRoot convention.
        reader >> *m12;  // Already swapped according to D* charge. m12 = m(pi+pi0)
        reader >> *m13;

        // Errors on Dalitz variables
        reader >> dummy;
        reader >> dummy;
        reader >> dummy;

        reader >> dummy; // Decay time
        reader >> dummy; // sigma_t

        reader >> dummy; // Md0
        reader >> dummy; // deltaM
        reader >> dummy; // ProbSig
        reader >> dummy; // Dst charge
        reader >> dummy; // Run
        reader >> dummy; // Event
        reader >> dummy; // Signal and four bkg fractions.
        reader >> dummy;
        reader >> dummy;
        reader >> dummy;
        reader >> dummy;

        // EXERCISE 1 (preliminary): Impose an artificial reconstruction efficiency
        // by throwing out events with a probability linear in m12.
        // NB! This exercise continues below.

        // EXERCISE 2: Instead of the above efficiency, impose a
        // K0 veto, by throwing out events with 0.475 < m23 < 0.505.

        // EXERCISE 3: Use both the above.

        eventNumber->setValue(data->getNumEvents());
        data->addEvent();

        dalitzplot.Fill(m12->getValue(), m13->getValue());
    }

    dalitzplot.SetStats(false);
    dalitzplot.Draw("colz");
    foodal->SaveAs("dalitzplot.png");
}

GooPdf *makeKzeroVeto() {
    if(kzero_veto)
        return kzero_veto;

    VetoInfo *kVetoInfo     = new VetoInfo();
    kVetoInfo->cyclic_index = PAIR_23;
    kVetoInfo->minimum      = new Variable("veto_min", 0.475 * 0.475);
    kVetoInfo->maximum      = new Variable("veto_max", 0.505 * 0.505);
    vector<VetoInfo *> vetos;
    vetos.push_back(kVetoInfo);
    kzero_veto = new DalitzVetoPdf("kzero_veto", m12, m13, motherM, neutrlM, chargeM, chargeM, vetos);
    return kzero_veto;
}

DalitzPlotPdf *makeSignalPdf(GooPdf *eff = 0) {
    DecayInfo *dtop0pp    = new DecayInfo();
    dtop0pp->motherMass   = _mD0;
    dtop0pp->daug1Mass    = piZeroMass;
    dtop0pp->daug2Mass    = piPlusMass;
    dtop0pp->daug3Mass    = piPlusMass;
    dtop0pp->meson_radius = 1.5;

    ResonancePdf *rhop = new ResonancePdf("rhop",
                                          ResPdfType::RBW,
                                          new Variable("rhop_amp_real", 1),
                                          new Variable("rhop_amp_imag", 0),
                                          fixedRhoMass,
                                          fixedRhoWidth,
                                          1,
                                          PAIR_12);

    bool fixAmps = false;

    ResonancePdf *rhom = new ResonancePdf(
        "rhom",
        ResPdfType::RBW,
        fixAmps ? new Variable("rhom_amp_real", 0.714) : new Variable("rhom_amp_real", 0.714, 0.001, 0, 0),
        fixAmps ? new Variable("rhom_amp_imag", -0.025) : new Variable("rhom_amp_imag", -0.025, 0.1, 0, 0),
        fixedRhoMass,
        fixedRhoWidth,
        1,
        PAIR_13);

    ResonancePdf *rho0 = new ResonancePdf(
        "rho0",
        ResPdfType::RBW,
        fixAmps ? new Variable("rho0_amp_real", 0.565) : new Variable("rho0_amp_real", 0.565, 0.001, 0, 0),
        fixAmps ? new Variable("rho0_amp_imag", 0.164) : new Variable("rho0_amp_imag", 0.164, 0.1, 0, 0),
        fixedRhoMass,
        fixedRhoWidth,
        1,
        PAIR_23);

    Variable *sharedMass = new Variable("rhop_1450_mass", 1.465, 0.01, 1.0, 2.0);
    Variable *shareWidth = new Variable("rhop_1450_width", 0.400, 0.01, 0.01, 5.0);

    ResonancePdf *rhop_1450 = new ResonancePdf(
        "rhop_1450",
        ResPdfType::RBW,
        fixAmps ? new Variable("rhop_1450_amp_real", -0.174) : new Variable("rhop_1450_amp_real", -0.174, 0.001, 0, 0),
        fixAmps ? new Variable("rhop_1450_amp_imag", -0.117) : new Variable("rhop_1450_amp_imag", -0.117, 0.1, 0, 0),
        sharedMass,
        shareWidth,
        1,
        PAIR_12);

    ResonancePdf *rho0_1450 = new ResonancePdf(
        "rho0_1450",
        ResPdfType::RBW,
        fixAmps ? new Variable("rho0_1450_amp_real", 0.325) : new Variable("rho0_1450_amp_real", 0.325, 0.001, 0, 0),
        fixAmps ? new Variable("rho0_1450_amp_imag", 0.057) : new Variable("rho0_1450_amp_imag", 0.057, 0.1, 0, 0),
        sharedMass,
        shareWidth,
        1,
        PAIR_23);

    ResonancePdf *rhom_1450 = new ResonancePdf(
        "rhom_1450",
        ResPdfType::RBW,
        fixAmps ? new Variable("rhom_1450_amp_real", 0.788) : new Variable("rhom_1450_amp_real", 0.788, 0.001, 0, 0),
        fixAmps ? new Variable("rhom_1450_amp_imag", 0.226) : new Variable("rhom_1450_amp_imag", 0.226, 0.1, 0, 0),
        sharedMass,
        shareWidth,
        1,
        PAIR_13);

    sharedMass = new Variable("rhop_1700_mass", 1.720, 0.01, 1.6, 1.9);
    shareWidth = new Variable("rhop_1700_width", 0.250, 0.01, 0.1, 1.0);

    ResonancePdf *rhop_1700 = new ResonancePdf(
        "rhop_1700",
        ResPdfType::RBW,
        fixAmps ? new Variable("rhop_1700_amp_real", 2.151) : new Variable("rhop_1700_amp_real", 2.151, 0.001, 0, 0),
        fixAmps ? new Variable("rhop_1700_amp_imag", -0.658) : new Variable("rhop_1700_amp_imag", -0.658, 0.1, 0, 0),
        sharedMass,
        shareWidth,
        1,
        PAIR_12);

    ResonancePdf *rho0_1700 = new ResonancePdf(
        "rho0_1700",
        ResPdfType::RBW,
        fixAmps ? new Variable("rho0_1700_amp_real", 2.400) : new Variable("rho0_1700_amp_real", 2.400, 0.001, 0, 0),
        fixAmps ? new Variable("rho0_1700_amp_imag", -0.734) : new Variable("rho0_1700_amp_imag", -0.734, 0.1, 0, 0),
        sharedMass,
        shareWidth,
        1,
        PAIR_23);

    ResonancePdf *rhom_1700 = new ResonancePdf(
        "rhom_1700",
        ResPdfType::RBW,
        fixAmps ? new Variable("rhom_1700_amp_real", 1.286) : new Variable("rhom_1700_amp_real", 1.286, 0.001, 0, 0),
        fixAmps ? new Variable("rhom_1700_amp_imag", -1.532) : new Variable("rhom_1700_amp_imag", -1.532, 0.1, 0, 0),
        sharedMass,
        shareWidth,
        1,
        PAIR_13);

    ResonancePdf *f0_980 = new ResonancePdf("f0_980",
                                            ResPdfType::RBW,
                                            fixAmps ? new Variable("f0_980_amp_real", 0.008 * (-_mD02))
                                                    : new Variable("f0_980_amp_real", 0.008 * (-_mD02), 0.001, 0, 0),
                                            fixAmps ? new Variable("f0_980_amp_imag", -0.013 * (-_mD02))
                                                    : new Variable("f0_980_amp_imag", -0.013 * (-_mD02), 0.1, 0, 0),
                                            new Variable("f0_980_mass", 0.980, 0.01, 0.8, 1.2),
                                            new Variable("f0_980_width", 0.044, 0.001, 0.001, 0.08),
                                            0,
                                            PAIR_23);

    ResonancePdf *f0_1370 = new ResonancePdf("f0_1370",
                                             ResPdfType::RBW,
                                             fixAmps ? new Variable("f0_1370_amp_real", -0.058 * (-_mD02))
                                                     : new Variable("f0_1370_amp_real", -0.058 * (-_mD02), 0.001, 0, 0),
                                             fixAmps ? new Variable("f0_1370_amp_imag", 0.026 * (-_mD02))
                                                     : new Variable("f0_1370_amp_imag", 0.026 * (-_mD02), 0.1, 0, 0),
                                             new Variable("f0_1370_mass", 1.434, 0.01, 1.2, 1.6),
                                             new Variable("f0_1370_width", 0.173, 0.01, 0.01, 0.4),
                                             0,
                                             PAIR_23);

    ResonancePdf *f0_1500 = new ResonancePdf("f0_1500",
                                             ResPdfType::RBW,
                                             fixAmps ? new Variable("f0_1500_amp_real", 0.057 * (-_mD02))
                                                     : new Variable("f0_1500_amp_real", 0.057 * (-_mD02), 0.001, 0, 0),
                                             fixAmps ? new Variable("f0_1500_amp_imag", 0.012 * (-_mD02))
                                                     : new Variable("f0_1500_amp_imag", 0.012 * (-_mD02), 0.1, 0, 0),
                                             new Variable("f0_1500_mass", 1.507, 0.01, 1.3, 1.7),
                                             new Variable("f0_1500_width", 0.109, 0.01, 0.01, 0.3),
                                             0,
                                             PAIR_23);

    ResonancePdf *f0_1710 = new ResonancePdf("f0_1710",
                                             ResPdfType::RBW,
                                             fixAmps ? new Variable("f0_1710_amp_real", 0.070 * (-_mD02))
                                                     : new Variable("f0_1710_amp_real", 0.070 * (-_mD02), 0.001, 0, 0),
                                             fixAmps ? new Variable("f0_1710_amp_imag", 0.087 * (-_mD02))
                                                     : new Variable("f0_1710_amp_imag", 0.087 * (-_mD02), 0.1, 0, 0),
                                             new Variable("f0_1710_mass", 1.714, 0.01, 1.5, 2.9),
                                             new Variable("f0_1710_width", 0.140, 0.01, 0.01, 0.5),
                                             0,
                                             PAIR_23);

    ResonancePdf *f2_1270
        = new ResonancePdf("f2_1270",
                           ResPdfType::RBW,
                           fixAmps ? new Variable("f2_1270_amp_real", -1.027 * (-_mD02inv))
                                   : new Variable("f2_1270_amp_real", -1.027 * (-_mD02inv), 0.001, 0, 0),
                           fixAmps ? new Variable("f2_1270_amp_imag", -0.162 * (-_mD02inv))
                                   : new Variable("f2_1270_amp_imag", -0.162 * (-_mD02inv), 0.1, 0, 0),
                           new Variable("f2_1270_mass", 1.2754, 0.01, 1.0, 1.5),
                           new Variable("f2_1270_width", 0.1851, 0.01, 0.01, 0.4),
                           2,
                           PAIR_23);

    ResonancePdf *f0_600 = new ResonancePdf("f0_600",
                                            ResPdfType::RBW,
                                            fixAmps ? new Variable("f0_600_amp_real", 0.068 * (-_mD02))
                                                    : new Variable("f0_600_amp_real", 0.068 * (-_mD02), 0.001, 0, 0),
                                            fixAmps ? new Variable("f0_600_amp_imag", 0.010 * (-_mD02))
                                                    : new Variable("f0_600_amp_imag", 0.010 * (-_mD02), 0.1, 0, 0),
                                            new Variable("f0_600_mass", 0.500, 0.01, 0.3, 0.7),
                                            new Variable("f0_600_width", 0.400, 0.01, 0.2, 0.6),
                                            0,
                                            PAIR_23);

    ResonancePdf *nonr = new ResonancePdf("nonr",
                                          ResPdfType::NONRES,
                                          fixAmps ? new Variable("nonr_amp_real", 0.5595 * (-1))
                                                  : new Variable("nonr_amp_real", 0.5595 * (-1), 0.001, 0, 0),
                                          fixAmps ? new Variable("nonr_amp_imag", -0.108761 * (-1))
                                                  : new Variable("nonr_amp_imag", -0.108761 * (-1), 0.1, 0, 0));

    dtop0pp->resonances.push_back(nonr);
    dtop0pp->resonances.push_back(rhop);
    dtop0pp->resonances.push_back(rho0);
    dtop0pp->resonances.push_back(rhom);
    dtop0pp->resonances.push_back(rhop_1450);
    dtop0pp->resonances.push_back(rho0_1450);
    dtop0pp->resonances.push_back(rhom_1450);
    dtop0pp->resonances.push_back(rhop_1700);
    dtop0pp->resonances.push_back(rho0_1700);
    dtop0pp->resonances.push_back(rhom_1700);
    dtop0pp->resonances.push_back(f0_980);
    dtop0pp->resonances.push_back(f0_1370);
    dtop0pp->resonances.push_back(f0_1500);
    dtop0pp->resonances.push_back(f0_1710);
    dtop0pp->resonances.push_back(f2_1270);
    dtop0pp->resonances.push_back(f0_600);

    if(!fitMasses) {
        for(vector<ResonancePdf *>::iterator res = dtop0pp->resonances.begin(); res != dtop0pp->resonances.end();
            ++res) {
            (*res)->setParameterConstantness(true);
        }
    }

    if(!eff) {
        // By default create a constant efficiency.
        vector<Variable *> offsets;
        vector<Variable *> observables;
        vector<Variable *> coefficients;

        observables.push_back(m12);
        observables.push_back(m13);
        offsets.push_back(constantZero);
        offsets.push_back(constantZero);
        coefficients.push_back(constantOne);
        eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0);
    }

    return new DalitzPlotPdf("signalPDF", m12, m13, eventNumber, dtop0pp, eff);
}

int runToyFit(std::string toyFileName, GooFit::Application &app) {
    m12 = new Variable("m12", 0, 3);
    m13 = new Variable("m13", 0, 3);
    m12->setNumBins(240);
    m13->setNumBins(240);
    eventNumber = new CountingVariable("eventNumber", 0, INT_MAX);
    getToyData(toyFileName, app);

    // EXERCISE 1 (real part): Create a PolynomialPdf which models
    // the efficiency you imposed in the preliminary, and use it in constructing
    // the signal PDF.

    // EXERCISE 2: Create a K0 veto function and use it as the efficiency.

    // EXERCISE 3: Make the efficiency a product of the two functions
    // from the previous exercises.

    DalitzPlotPdf *signal = makeSignalPdf();
    signal->setData(data);
    signal->setDataSize(data->getNumEvents());
    FitManager datapdf(signal);

    datapdf.fit();
    return datapdf;
}

int main(int argc, char **argv) {
    GooFit::Application app("Dalitz example", argc, argv);

    std::string filename = "dalitz_toyMC_000.txt";
    app.add_option("-f,--filename,filename", filename, "File to read in", true)->check(GooFit::ExistingFile);

    try {
        app.run();
    } catch(const GooFit::ParseError &e) {
        return app.exit(e);
    }

    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasColor(10);
    gStyle->SetFrameFillColor(10);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetPadColor(0);
    gStyle->SetTitleColor(1);
    gStyle->SetStatColor(0);
    gStyle->SetFillColor(0);
    gStyle->SetFuncWidth(1);
    gStyle->SetLineWidth(1);
    gStyle->SetLineColor(1);
    gStyle->SetPalette(1, 0);
    foo    = new TCanvas();
    foodal = new TCanvas();
    foodal->Size(10, 10);

    try {
        return runToyFit(filename, app);
    } catch(const std::runtime_error &e) {
        std::cerr << e.what() << std::endl;
        return 7;
    }
}
