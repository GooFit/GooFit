// ROOT stuff
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TRandom.h>
#include <TRandom3.h>

// System stuff
#include <fstream>
#include <sys/time.h>
#include <sys/times.h>

// GooFit stuff
#include <goofit/Application.h>
#include <goofit/FitManager.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/basic/PolynomialPdf.h>
#include <goofit/PDFs/combine/AddPdf.h>
#include <goofit/PDFs/combine/ProdPdf.h>
#include <goofit/PDFs/physics/DalitzPlotPdf.h>
#include <goofit/PDFs/physics/DalitzPlotter.h>
#include <goofit/PDFs/physics/DalitzVetoPdf.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>
#include <goofit/detail/Style.h>

using namespace std;
using namespace GooFit;

Variable fixedRhoMass("rho_mass", 0.7758, 0.01, 0.7, 0.8);
Variable fixedRhoWidth("rho_width", 0.1503, 0.01, 0.1, 0.2);

const fptype _mD0       = 1.86484;
const fptype _mD02      = _mD0 * _mD0;
const fptype _mD02inv   = 1. / _mD02;
const fptype piPlusMass = 0.13957018;
const fptype piZeroMass = 0.1349766;

// Constants used in more than one PDF component.
Variable motherM("motherM", _mD0);
Variable chargeM("chargeM", piPlusMass);
Variable neutrlM("neutrlM", piZeroMass);
Variable massSum("massSum", _mD0 *_mD0 + 2 * piPlusMass * piPlusMass + piZeroMass * piZeroMass); // = 3.53481
Variable constantOne("constantOne", 1);
Variable constantZero("constantZero", 0);

fptype cpuGetM23(fptype massPZ, fptype massPM) {
    return (_mD02 + piZeroMass * piZeroMass + piPlusMass * piPlusMass + piPlusMass * piPlusMass - massPZ - massPM);
}

void getToyData(std::string toyFileName, GooFit::Application &app, DataSet &data) {
    toyFileName = app.get_filename(toyFileName, "examples/dalitz");

    auto obs               = data.getObservables();
    Observable m12         = obs.at(0);
    Observable m13         = obs.at(1);
    Observable eventNumber = obs.at(2);

    TH2F dalitzplot("dalitzplot",
                    "Original Data",
                    m12.getNumBins(),
                    m12.getLowerLimit(),
                    m12.getUpperLimit(),
                    m13.getNumBins(),
                    m13.getLowerLimit(),
                    m13.getUpperLimit());
    std::vector<Observable> vars;
    vars.push_back(m12);
    vars.push_back(m13);
    vars.push_back(eventNumber);

    std::ifstream reader(toyFileName);
    std::string buffer;

    while(!reader.eof()) {
        reader >> buffer;

        if(buffer == "====")
            break;

        std::cout << buffer;
    }

    double dummy = 0;

    while(reader >> dummy
          // m23, m(pi+ pi-), called m12 in processToyRoot convention.
          // Already swapped according to D* charge. m12 = m(pi+pi0)
          >> dummy >> m12 >> m13

          // Errors on Dalitz variables
          >> dummy >> dummy >> dummy

          // Decay time, sigma_t
          >> dummy >> dummy

          // Md0, deltaM, Prob, Sig, Dst charge, Run, Event, Signal and four bkg fractions
          >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy) {
        // EXERCISE 1 (preliminary): Impose an artificial reconstruction efficiency
        // by throwing out events with a probability linear in m12.
        // NB! This exercise continues below.

        // EXERCISE 2: Instead of the above efficiency, impose a
        // K0 veto, by throwing out events with 0.475 < m23 < 0.505.

        // EXERCISE 3: Use both the above.

        eventNumber.setValue(data.getNumEvents());
        data.addEvent();

        dalitzplot.Fill(m12.getValue(), m13.getValue());
    }

    GOOFIT_INFO("Read in {} events", data.getNumEvents());

    TCanvas foo;
    dalitzplot.SetStats(false);
    dalitzplot.Draw("colz");
    foo.SaveAs("dalitzplot.png");
}

void makeToyData(DalitzPlotter &dplotter, UnbinnedDataSet &data) {}

DalitzPlotPdf *makeSignalPdf(Observable m12, Observable m13, EventNumber eventNumber, GooPdf *eff = 0) {
    DecayInfo3 dtop0pp;
    dtop0pp.motherMass   = _mD0;
    dtop0pp.daug1Mass    = piZeroMass;
    dtop0pp.daug2Mass    = piPlusMass;
    dtop0pp.daug3Mass    = piPlusMass;
    dtop0pp.meson_radius = 1.5;

    ResonancePdf *rhop = new Resonances::RBW(
        "rhop", Variable("rhop_amp_real", 1), Variable("rhop_amp_imag", 0), fixedRhoMass, fixedRhoWidth, 1, PAIR_12);

    bool fixAmps = false; // Takes ~400x longer

    ResonancePdf *rhom = new Resonances::RBW(
        "rhom",
        fixAmps ? Variable("rhom_amp_real", 0.714) : Variable("rhom_amp_real", 0.714, 0.001, 0, 0),
        fixAmps ? Variable("rhom_amp_imag", -0.025) : Variable("rhom_amp_imag", -0.025, 0.1, 0, 0),
        fixedRhoMass,
        fixedRhoWidth,
        1,
        PAIR_13);

    ResonancePdf *rho0 = new Resonances::RBW(
        "rho0",
        fixAmps ? Variable("rho0_amp_real", 0.565) : Variable("rho0_amp_real", 0.565, 0.001, 0, 0),
        fixAmps ? Variable("rho0_amp_imag", 0.164) : Variable("rho0_amp_imag", 0.164, 0.1, 0, 0),
        fixedRhoMass,
        fixedRhoWidth,
        1,
        PAIR_23);

    Variable sharedMass("rhop_1450_mass", 1.465, 0.01, 1.0, 2.0);
    Variable shareWidth("rhop_1450_width", 0.400, 0.01, 0.01, 5.0);

    ResonancePdf *rhop_1450 = new Resonances::RBW(
        "rhop_1450",
        fixAmps ? Variable("rhop_1450_amp_real", -0.174) : Variable("rhop_1450_amp_real", -0.174, 0.001, 0, 0),
        fixAmps ? Variable("rhop_1450_amp_imag", -0.117) : Variable("rhop_1450_amp_imag", -0.117, 0.1, 0, 0),
        sharedMass,
        shareWidth,
        1,
        PAIR_12);

    ResonancePdf *rho0_1450 = new Resonances::RBW(
        "rho0_1450",
        fixAmps ? Variable("rho0_1450_amp_real", 0.325) : Variable("rho0_1450_amp_real", 0.325, 0.001, 0, 0),
        fixAmps ? Variable("rho0_1450_amp_imag", 0.057) : Variable("rho0_1450_amp_imag", 0.057, 0.1, 0, 0),
        sharedMass,
        shareWidth,
        1,
        PAIR_23);

    ResonancePdf *rhom_1450 = new Resonances::RBW(
        "rhom_1450",
        fixAmps ? Variable("rhom_1450_amp_real", 0.788) : Variable("rhom_1450_amp_real", 0.788, 0.001, 0, 0),
        fixAmps ? Variable("rhom_1450_amp_imag", 0.226) : Variable("rhom_1450_amp_imag", 0.226, 0.1, 0, 0),
        sharedMass,
        shareWidth,
        1,
        PAIR_13);

    Variable sharedMass2("rhop_1700_mass", 1.720, 0.01, 1.6, 1.9);
    Variable shareWidth2("rhop_1700_width", 0.250, 0.01, 0.1, 1.0);

    ResonancePdf *rhop_1700 = new Resonances::RBW(
        "rhop_1700",
        fixAmps ? Variable("rhop_1700_amp_real", 2.151) : Variable("rhop_1700_amp_real", 2.151, 0.001, 0, 0),
        fixAmps ? Variable("rhop_1700_amp_imag", -0.658) : Variable("rhop_1700_amp_imag", -0.658, 0.1, 0, 0),
        sharedMass2,
        shareWidth2,
        1,
        PAIR_12);

    ResonancePdf *rho0_1700 = new Resonances::RBW(
        "rho0_1700",
        fixAmps ? Variable("rho0_1700_amp_real", 2.400) : Variable("rho0_1700_amp_real", 2.400, 0.001, 0, 0),
        fixAmps ? Variable("rho0_1700_amp_imag", -0.734) : Variable("rho0_1700_amp_imag", -0.734, 0.1, 0, 0),
        sharedMass2,
        shareWidth2,
        1,
        PAIR_23);

    ResonancePdf *rhom_1700 = new Resonances::RBW(
        "rhom_1700",
        fixAmps ? Variable("rhom_1700_amp_real", 1.286) : Variable("rhom_1700_amp_real", 1.286, 0.001, 0, 0),
        fixAmps ? Variable("rhom_1700_amp_imag", -1.532) : Variable("rhom_1700_amp_imag", -1.532, 0.1, 0, 0),
        sharedMass2,
        shareWidth2,
        1,
        PAIR_13);

    ResonancePdf *f0_980 = new Resonances::RBW("f0_980",
                                               fixAmps ? Variable("f0_980_amp_real", 0.008 * (-_mD02))
                                                       : Variable("f0_980_amp_real", 0.008 * (-_mD02), 0.001, 0, 0),
                                               fixAmps ? Variable("f0_980_amp_imag", -0.013 * (-_mD02))
                                                       : Variable("f0_980_amp_imag", -0.013 * (-_mD02), 0.1, 0, 0),
                                               Variable("f0_980_mass", 0.980, 0.01, 0.8, 1.2),
                                               Variable("f0_980_width", 0.044, 0.001, 0.001, 0.08),
                                               0,
                                               PAIR_23);

    ResonancePdf *f0_1370 = new Resonances::RBW("f0_1370",
                                                fixAmps ? Variable("f0_1370_amp_real", -0.058 * (-_mD02))
                                                        : Variable("f0_1370_amp_real", -0.058 * (-_mD02), 0.001, 0, 0),
                                                fixAmps ? Variable("f0_1370_amp_imag", 0.026 * (-_mD02))
                                                        : Variable("f0_1370_amp_imag", 0.026 * (-_mD02), 0.1, 0, 0),
                                                Variable("f0_1370_mass", 1.434, 0.01, 1.2, 1.6),
                                                Variable("f0_1370_width", 0.173, 0.01, 0.01, 0.4),
                                                0,
                                                PAIR_23);

    ResonancePdf *f0_1500 = new Resonances::RBW("f0_1500",
                                                fixAmps ? Variable("f0_1500_amp_real", 0.057 * (-_mD02))
                                                        : Variable("f0_1500_amp_real", 0.057 * (-_mD02), 0.001, 0, 0),
                                                fixAmps ? Variable("f0_1500_amp_imag", 0.012 * (-_mD02))
                                                        : Variable("f0_1500_amp_imag", 0.012 * (-_mD02), 0.1, 0, 0),
                                                Variable("f0_1500_mass", 1.507, 0.01, 1.3, 1.7),
                                                Variable("f0_1500_width", 0.109, 0.01, 0.01, 0.3),
                                                0,
                                                PAIR_23);

    ResonancePdf *f0_1710 = new Resonances::RBW("f0_1710",
                                                fixAmps ? Variable("f0_1710_amp_real", 0.070 * (-_mD02))
                                                        : Variable("f0_1710_amp_real", 0.070 * (-_mD02), 0.001, 0, 0),
                                                fixAmps ? Variable("f0_1710_amp_imag", 0.087 * (-_mD02))
                                                        : Variable("f0_1710_amp_imag", 0.087 * (-_mD02), 0.1, 0, 0),
                                                Variable("f0_1710_mass", 1.714, 0.01, 1.5, 2.9),
                                                Variable("f0_1710_width", 0.140, 0.01, 0.01, 0.5),
                                                0,
                                                PAIR_23);

    ResonancePdf *f2_1270
        = new Resonances::RBW("f2_1270",
                              fixAmps ? Variable("f2_1270_amp_real", -1.027 * (-_mD02inv))
                                      : Variable("f2_1270_amp_real", -1.027 * (-_mD02inv), 0.001, 0, 0),
                              fixAmps ? Variable("f2_1270_amp_imag", -0.162 * (-_mD02inv))
                                      : Variable("f2_1270_amp_imag", -0.162 * (-_mD02inv), 0.1, 0, 0),
                              Variable("f2_1270_mass", 1.2754, 0.01, 1.0, 1.5),
                              Variable("f2_1270_width", 0.1851, 0.01, 0.01, 0.4),
                              2,
                              PAIR_23);

    ResonancePdf *f0_600 = new Resonances::RBW("f0_600",
                                               fixAmps ? Variable("f0_600_amp_real", 0.068 * (-_mD02))
                                                       : Variable("f0_600_amp_real", 0.068 * (-_mD02), 0.001, 0, 0),
                                               fixAmps ? Variable("f0_600_amp_imag", 0.010 * (-_mD02))
                                                       : Variable("f0_600_amp_imag", 0.010 * (-_mD02), 0.1, 0, 0),
                                               Variable("f0_600_mass", 0.500, 0.01, 0.3, 0.7),
                                               Variable("f0_600_width", 0.400, 0.01, 0.2, 0.6),
                                               0,
                                               PAIR_23);

    ResonancePdf *nonr = new Resonances::NonRes(
        "nonr",
        fixAmps ? Variable("nonr_amp_real", 0.5595 * (-1)) : Variable("nonr_amp_real", 0.5595 * (-1), 0.001, 0, 0),
        fixAmps ? Variable("nonr_amp_imag", -0.108761 * (-1)) : Variable("nonr_amp_imag", -0.108761 * (-1), 0.1, 0, 0));

    dtop0pp.resonances.push_back(nonr);
    dtop0pp.resonances.push_back(rhop);
    dtop0pp.resonances.push_back(rho0);
    dtop0pp.resonances.push_back(rhom);
    dtop0pp.resonances.push_back(rhop_1450);
    dtop0pp.resonances.push_back(rho0_1450);
    dtop0pp.resonances.push_back(rhom_1450);
    dtop0pp.resonances.push_back(rhop_1700);
    dtop0pp.resonances.push_back(rho0_1700);
    dtop0pp.resonances.push_back(rhom_1700);
    dtop0pp.resonances.push_back(f0_980);
    dtop0pp.resonances.push_back(f0_1370);
    dtop0pp.resonances.push_back(f0_1500);
    dtop0pp.resonances.push_back(f0_1710);
    dtop0pp.resonances.push_back(f2_1270);
    dtop0pp.resonances.push_back(f0_600);

    bool fitMasses = false;

    if(!fitMasses) {
        for(vector<ResonancePdf *>::iterator res = dtop0pp.resonances.begin(); res != dtop0pp.resonances.end(); ++res) {
            (*res)->setParameterConstantness(true);
        }
    }

    if(!eff) {
        // By default create a constant efficiency.
        vector<Variable> offsets       = {constantZero, constantZero};
        vector<Observable> observables = {m12, m13};
        vector<Variable> coefficients  = {constantOne};

        eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0);
    }

    return new DalitzPlotPdf("signalPDF", m12, m13, eventNumber, dtop0pp, eff);
}

int runToyFit(DalitzPlotPdf *signal, UnbinnedDataSet *data) {
    // EXERCISE 1 (real part): Create a PolynomialPdf which models
    // the efficiency you imposed in the preliminary, and use it in constructing
    // the signal PDF.

    // EXERCISE 2: Create a K0 veto function and use it as the efficiency.

    // EXERCISE 3: Make the efficiency a product of the two functions
    // from the previous exercises.

    signal->setData(data);
    signal->setDataSize(data->getNumEvents());
    FitManager datapdf(signal);

    datapdf.fit();

    ProdPdf prodpdf{"prodpdf", {signal}};

    DalitzPlotter plotter(&prodpdf, signal);

    TCanvas foo;
    TH2F *dalitzplot = plotter.make2D();
    dalitzplot->Draw("colz");

    foo.SaveAs("dalitzpdf.png");

    return datapdf;
}

int main(int argc, char **argv) {
    GooFit::Application app("Dalitz example", argc, argv);

    std::string filename = "dalitz_toyMC_000.txt";
    app.add_option("-f,--filename,filename", filename, "File to read in", true)->check(GooFit::ExistingFile);

    bool make_toy;
    app.add_flag("-m,--make-toy", make_toy, "Make a toy instead of reading a file in");

    GOOFIT_PARSE(app);

    GooFit::setROOTStyle();

    // Observables setup
    Observable m12("m12", 0, 3);
    Observable m13("m13", 0, 3);
    EventNumber eventNumber("eventNumber");
    m12.setNumBins(240);
    m13.setNumBins(240);

    // Prepare the data
    UnbinnedDataSet data({m12, m13, eventNumber});

    // Set up the model
    DalitzPlotPdf *signal = makeSignalPdf(m12, m13, eventNumber);

    // A wrapper for plotting without complex number segfault
    ProdPdf prodpdf{"prodpdf", {signal}};

    // Add nice tool for making data or plotting
    DalitzPlotter dplotter{&prodpdf, signal};

    // Read in data
    if(make_toy) {
        dplotter.fillDataSetMC(data, 1000000);
    } else {
        getToyData(filename, app, data);
    }

    try {
        return runToyFit(signal, &data);
    } catch(const std::runtime_error &e) {
        std::cerr << e.what() << std::endl;
        return 7;
    }
}
