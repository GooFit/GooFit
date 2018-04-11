// ROOT stuff
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TText.h>

// System stuff
#include <cassert>
#include <climits>
#include <fstream>

// GooFit stuff
#include <goofit/Application.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/basic/BinTransformPdf.h>
#include <goofit/PDFs/basic/CrystalBallPdf.h>
#include <goofit/PDFs/basic/ExpGausPdf.h>
#include <goofit/PDFs/basic/ExpPdf.h>
#include <goofit/PDFs/basic/GaussianPdf.h>
#include <goofit/PDFs/basic/JohnsonSUPdf.h>
#include <goofit/PDFs/basic/PolynomialPdf.h>
#include <goofit/PDFs/basic/SmoothHistogramPdf.h>
#include <goofit/PDFs/basic/StepPdf.h>
#include <goofit/PDFs/basic/TrigThresholdPdf.h>
#include <goofit/PDFs/basic/VoigtianPdf.h>
#include <goofit/PDFs/combine/AddPdf.h>
#include <goofit/PDFs/combine/ConvolutionPdf.h>
#include <goofit/PDFs/combine/EventWeightedAddPdf.h>
#include <goofit/PDFs/combine/MappedPdf.h>
#include <goofit/PDFs/combine/ProdPdf.h>
#include <goofit/PDFs/physics/DalitzVetoPdf.h>
#include <goofit/PDFs/physics/IncoherentSumPdf.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/PDFs/physics/TddpPdf.h>
#include <goofit/PDFs/physics/ThreeGaussResolution_Aux.h>
#include <goofit/PDFs/physics/TruthResolution_Aux.h>
#include <goofit/Variable.h>
#include <goofit/detail/Uncertain.h>

#include <goofit/fitting/FitManagerMinuit1.h>
#include <goofit/fitting/FitManagerMinuit2.h>

#include <goofit/FunctorWriter.h>
#include <goofit/PDFs/combine/CompositePdf.h>
#include <goofit/UnbinnedDataSet.h>

using namespace std;
using namespace GooFit;

TCanvas *foo;
TCanvas *foodal;

UnbinnedDataSet *data     = nullptr;
UnbinnedDataSet *effdata  = nullptr;
BinnedDataSet *binEffData = nullptr;
TH2F *weightHistogram     = nullptr;
TH2F *underlyingBins      = nullptr;

// I hate having to use globals, but this is the best way for now
GooFit::Application *app_ptr;
bool minuit1;

Observable *m12          = nullptr;
Observable *m13          = nullptr;
EventNumber *eventNumber = nullptr;
Observable *massd0       = nullptr;
Observable *deltam       = nullptr;
Observable *dtime        = nullptr;
Observable *sigma        = nullptr;
Observable *wSig0        = nullptr;
Observable *wBkg1        = nullptr;
Observable *wBkg2        = nullptr;
Observable *wBkg3        = nullptr;
Observable *wBkg4        = nullptr;

bool fitMasses = false;
Variable fixedRhoMass("rho_mass", 0.7758, 0.01, 0.7, 0.8);
Variable fixedRhoWidth("rho_width", 0.1503, 0.01, 0.1, 0.2);

// Systematic variables
double luckyFrac           = 0.5;
double mesonRad            = 1.5;
int normBinning            = 240;
int blindSeed              = 4;
int mdslices               = 1;
double md0offset           = 0;
int md0_lower_window       = -2;
int md0_upper_window       = 2;
double deltam_lower_window = -2;
double deltam_upper_window = 2;
double lowerTime           = -2;
double upperTime           = 3;
double maxSigma            = 0.8;
bool polyEff               = false;
bool saveEffPlot           = true;
bool useHistogramSigma     = false;
enum Bkg2Model { Histogram, Parameter, Sideband };
Bkg2Model bkg2Model        = Sideband;
std::string bkg2Model_str  = "sideband";
bool notUseBackground3Hist = false;
bool notUseBackground4Hist = false;
bool makePlots             = false;
int m23Slices              = 6;
bool drop_rho_1450         = false;
bool drop_rho_1700         = false;
bool drop_f0_980           = false;
bool drop_f0_1370          = false;
bool drop_f0_1500          = false;
bool drop_f0_1710          = false;
bool drop_f2_1270          = false;
bool drop_f0_600           = false;
bool gaussBkgTime          = false;
bool mikhailSetup          = false;
int bkgHistBins            = 80;
string paramUp             = "";
string paramDn             = "";
int bkgHistRandSeed        = -1;

const fptype _mD0       = 1.86484;
const fptype _mD02      = _mD0 * _mD0;
const fptype _mD02inv   = 1. / _mD02;
const fptype piPlusMass = 0.13957018;
const fptype piZeroMass = 0.1349766;
enum Resolutions { DplotRes = 1, TimeRes = 2, Efficiency = 4 };

size_t maxEvents = 100000;

Variable motherM("motherM", 1.86484);
Variable chargeM("chargeM", 0.13957018);
Variable neutrlM("neutrlM", 0.1349766);
Variable constantBigM("constantBigM", _mD02 + 2 * piPlusMass * piPlusMass + piZeroMass * piZeroMass);

// Constants used in more than one PDF component.
Variable constantOne("constantOne", 1);
Variable constantTwo("constantTwo", 2);
Variable constantZero("constantZero", 0);
Variable constantMinusOne("constantMinusOne", -1);
Variable minDalitzX("minDalitzX", pow(piPlusMass + piZeroMass, 2));
Variable maxDalitzX("maxDalitzX", pow(_mD0 - piPlusMass, 2));
Variable minDalitzY("minDalitzY", pow(piPlusMass + piZeroMass, 2));
Variable maxDalitzY("maxDalitzY", pow(_mD0 - piPlusMass, 2));
Variable minDalitzZ("minDalitzZ", pow(piPlusMass + piPlusMass, 2));
Variable maxDalitzZ("maxDalitzZ", pow(_mD0 - piZeroMass, 2));

std::vector<Variable> weights;
std::vector<Observable> obsweights;

std::vector<PdfBase *> comps;
TH1F *dataTimePlot        = nullptr;
TH1F *loM23Sigma          = nullptr;
TH1F *hiM23Sigma          = nullptr;
TddpPdf *signalDalitz     = nullptr;
IncoherentSumPdf *incsum1 = nullptr;
IncoherentSumPdf *incsum2 = nullptr;
IncoherentSumPdf *incsum3 = nullptr;
IncoherentSumPdf *incsum4 = nullptr;
IncoherentSumPdf *incsum5 = nullptr;
IncoherentSumPdf *incsum6 = nullptr;
GooPdf *sig0_jsugg        = nullptr;
GooPdf *bkg2_jsugg        = nullptr;
GooPdf *bkg3_jsugg        = nullptr;
GooPdf *bkg4_jsugg        = nullptr;
Variable massSum("massSum", _mD0 *_mD0 + 2 * piPlusMass * piPlusMass + piZeroMass * piZeroMass); // = 3.53481
GooPdf *kzero_veto = nullptr;

bool doToyStudy            = false;
float md0_toy_width        = 0.0075;
float md0_toy_mean         = 1.8654;
const float toySigFraction = 0.6;
const float toyBkgTimeMean = 0.0;
const float toyBkgTimeRMS  = 0.7;
std::string toyFileName;
char strbuffer[1000];

SmoothHistogramPdf *makeBackgroundHistogram(int bkgnum, std::string overridename = "");
void makeToyDalitzPlots(GooPdf *overallSignal, std::string plotdir = "./plots_from_toy_mixfit/");
void getBackgroundFile(int bkgType);

double intGaus = -1;
double calcGauInteg(double x1, double x2) {
    if(x1 > x2) {
        double swp = x2;
        x2         = x1;
        x1         = swp;
    }

    double sum1 = x1, sum2 = x2;
    double value1 = x1, value2 = x2;

    for(int i = 1; i < 100; i++) {
        value1 = value1 * x1 * x1 / (2 * i + 1);
        sum1 += value1;
        value2 = value2 * x2 * x2 / (2 * i + 1);
        sum2 += value2;
    }

    return sum2 * exp(-(x2 * x2) / 2) - sum1 * exp(-(x1 * x1) / 2);
}

double calcToyWeight(double sigratio, double m) {
    if(intGaus < 0)
        intGaus = calcGauInteg(0.0075 * md0_lower_window / md0_toy_width, 0.0075 * md0_upper_window / md0_toy_width);

    double t    = (m - md0_toy_mean) / md0_toy_width;
    double fsig = sigratio * exp(-t * t / 2.) / intGaus;
    double fbkg = (1 - sigratio) / ((md0_upper_window - md0_lower_window) * 0.0075 / md0_toy_width);
    return fsig / (fsig + fbkg);
}

void printMemoryStatus(std::string file, int line);
void loadDataFile(std::string fname, UnbinnedDataSet **setToFill = 0, int effSkip = 3);
int runBackgroundDalitzFit(int bkgType, bool plots = false);

void normalize(TH1F *dat) {
    double integral = 0;

    for(int i = 1; i <= dat->GetNbinsX(); ++i) {
        integral += dat->GetBinContent(i);
    }

    integral = 1.0 / integral;

    for(int i = 1; i <= dat->GetNbinsX(); ++i) {
        dat->SetBinContent(i, integral * dat->GetBinContent(i));
        dat->SetBinError(i, integral * dat->GetBinError(i));
    }
}

fptype cpuGetM23(fptype massPZ, fptype massPM) {
    return (_mD02 + piZeroMass * piZeroMass + piPlusMass * piPlusMass + piPlusMass * piPlusMass - massPZ - massPM);
}

bool cpuDalitz(fptype m12, fptype m13, fptype bigM, fptype dm1, fptype dm2, fptype dm3) {
    if(m12 < pow(dm1 + dm2, 2))
        return false; // This m12 cannot exist, it's less than the square of the (1,2) particle mass.

    if(m12 > pow(bigM - dm3, 2))
        return false; // This doesn't work either, there's no room for an at-rest 3 daughter.

    // Calculate energies of 1 and 3 particles in m12 rest frame.
    fptype e1star = 0.5 * (m12 - dm2 * dm2 + dm1 * dm1) / sqrt(m12);
    fptype e3star = 0.5 * (bigM * bigM - m12 - dm3 * dm3) / sqrt(m12);

    // Bounds for m13 at this value of m12.
    fptype minimum
        = pow(e1star + e3star, 2) - pow(sqrt(e1star * e1star - dm1 * dm1) + sqrt(e3star * e3star - dm3 * dm3), 2);

    if(m13 < minimum)
        return false;

    fptype maximum
        = pow(e1star + e3star, 2) - pow(sqrt(e1star * e1star - dm1 * dm1) - sqrt(e3star * e3star - dm3 * dm3), 2);

    if(m13 > maximum)
        return false;

    return true;
}

void plotLoHiSigma() {
    if(!loM23Sigma)
        return;

    normalize(loM23Sigma);
    normalize(hiM23Sigma);

    loM23Sigma->SetLineWidth(3);
    loM23Sigma->SetLineColor(kBlue);
    hiM23Sigma->SetLineWidth(3);
    hiM23Sigma->SetLineColor(kRed);

    hiM23Sigma->GetXaxis()->SetTitle("#sigma_{t} [ps]");
    hiM23Sigma->GetYaxis()->SetTitle("Fraction / 8 fs");
    hiM23Sigma->Draw("l");
    loM23Sigma->Draw("lsame");

    TLegend blah(0.7, 0.8, 0.95, 0.95);
    blah.AddEntry(loM23Sigma, "m_{0}^{2} < 1.5", "l");
    blah.AddEntry(hiM23Sigma, "m_{0}^{2} > 1.5", "l");
    blah.Draw();

    foo->SaveAs("./plots_from_mixfit/sigmacomp.png");
}

void plotFit(Observable var, UnbinnedDataSet *dat, GooPdf *fit) {
    int numEvents  = dat->getNumEvents();
    TH1F *dat_hist = new TH1F(
        (var.getName() + "_dathist").c_str(), "", var.getNumBins(), var.getLowerLimit(), var.getUpperLimit());

    for(int i = 0; i < numEvents; ++i) {
        dat_hist->Fill(dat->getValue(var, i));
    }

    TH1F *pdf_hist = new TH1F(
        (var.getName() + "_pdfhist").c_str(), "", var.getNumBins(), var.getLowerLimit(), var.getUpperLimit());
    std::vector<fptype> values = fit->evaluateAtPoints(var);

    double totalPdf = 0;

    for(int i = 0; i < var.getNumBins(); ++i) {
        totalPdf += values[i];
    }

    for(int i = 0; i < var.getNumBins(); ++i) {
        pdf_hist->SetBinContent(i + 1, values[i] * numEvents / totalPdf);
    }

    pdf_hist->SetStats(false);
    pdf_hist->SetLineColor(kBlue);
    pdf_hist->SetLineWidth(3);

    dat_hist->SetStats(false);
    dat_hist->SetMarkerStyle(8);
    dat_hist->SetMarkerSize(0.7);
    dat_hist->Draw("p");

    pdf_hist->Draw("lsame");

    foo->SaveAs((var.getName() + "_fit.png").c_str());
    foo->SetLogy(true);
    foo->SaveAs((var.getName() + "_logfit.png").c_str());
    foo->SetLogy(false);
}

bool readWrapper(std::ifstream &reader, std::string fname = strbuffer) {
    std::cout << "Now open file " << fname << " for reading" << std::endl;
    reader.open(fname.c_str());
    assert(reader.good());
    return true;
}

void getToyData(float sigweight = 0.9) {
    if(!data) {
        std::vector<Observable> vars;
        vars.push_back(*m12);
        vars.push_back(*m13);
        vars.push_back(*dtime);
        vars.push_back(*sigma);
        vars.push_back(*eventNumber);
        vars.push_back(*wSig0);
        //  vars.push_back(wBkg1);
        //  vars.push_back(wBkg2);
        data = new UnbinnedDataSet(vars);
    }

    std::ifstream reader;
    readWrapper(reader, toyFileName);
    std::string buffer;

    while(!reader.eof()) {
        reader >> buffer;

        if(buffer == "====")
            break;

        std::cout << buffer;
    }

    TRandom3 donram(42);

    int nsig       = 0;
    double sigprob = 0;
    double dummy   = 0;
    double md0     = md0_toy_mean;

    while(reader >> dummy

          // Ignoring m23, m(pi+ pi-), called m12 in processToyRoot convention.
          // Already swapped m12 m13 according to D* charge. m12 = m(pi+pi0)
          >> dummy >> *m12 >> *m13

          // Errors on Dalitz variables
          >> dummy >> dummy >> dummy

          // Remaining two interesting values
          >> *dtime >> *sigma

          // Md0, deltaM, ProbSig, Dst charge, Run, Event, Signal and four bkg fractions
          >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy) {
        double resolution = donram.Gaus(0, 1);
        dtime->setValue(dtime->getValue() + resolution * sigma->getValue());

        eventNumber->setValue(data->getNumEvents());

        do {
            md0 = donram.Gaus(md0_toy_mean, md0_toy_width);
        } while(md0 <= 1.8654 + 0.0075 * md0_lower_window || md0 >= 1.8654 + 0.0075 * md0_upper_window);

        //    wSig0->getValue() = sigweight;
        wSig0->setValue(calcToyWeight(sigweight, md0));
        sigprob += wSig0->getValue();
        data->addEvent();
        nsig++;

        if(data->getNumEvents() < 10) {
            std::cout << data->getNumEvents() << " : " << m12->getValue() << " " << m13->getValue() << " "
                      << dtime->getValue() << " " << sigma->getValue() << " " << std::endl;
        }

        if(data->getNumEvents() >= maxEvents)
            break;
    }

    GOOFIT_INFO("Read in {} events", data->getNumEvents())

    for(int ib = 0; ib < nsig * (1 - sigweight) / sigweight; ib++) {
        do {
            m12->setValue(donram.Uniform() * (m12->getUpperLimit() - m12->getLowerLimit()) + m12->getLowerLimit());
            m13->setValue(donram.Uniform() * (m13->getUpperLimit() - m13->getLowerLimit()) + m13->getLowerLimit());
        } while(!cpuDalitz(m12->getValue(), m13->getValue(), _mD0, piZeroMass, piPlusMass, piPlusMass));

        do {
            dtime->setValue(donram.Gaus(toyBkgTimeMean, toyBkgTimeRMS));
        } while(!(dtime->getValue() > dtime->getLowerLimit() && dtime->getValue() < dtime->getUpperLimit()));

        eventNumber->setValue(data->getNumEvents());
        md0 = donram.Uniform(1.8654 + 0.0075 * md0_lower_window, 1.8654 + 0.0075 * md0_upper_window);
        //    wSig0->getValue() = sigweight;
        wSig0->setValue(calcToyWeight(sigweight, md0));
        sigprob += wSig0->getValue();
        data->addEvent();
    }

    reader.close();
}

GooPdf *makeEfficiencyPdf() {
    Variable effSmoothing("effSmoothing", 1.0, 0.1, 0, 1.25);
    // Variable* effSmoothing = new Variable("effSmoothing", 0);
    SmoothHistogramPdf *ret = new SmoothHistogramPdf("efficiency", binEffData, effSmoothing);
    return ret;
}

GooPdf *makeKzeroVeto() {
    if(kzero_veto)
        return kzero_veto;

    VetoInfo kVetoInfo(Variable("veto_min", 0.475 * 0.475), Variable("veto_max", 0.505 * 0.505), PAIR_23);

    vector<VetoInfo> vetos;
    vetos.push_back(kVetoInfo);
    kzero_veto = new DalitzVetoPdf("kzero_veto", *m12, *m13, motherM, neutrlM, chargeM, chargeM, vetos);
    return kzero_veto;
}

GooPdf *makeEfficiencyPoly() {
    if(!kzero_veto)
        makeKzeroVeto();

    vector<Variable> offsets;
    vector<Observable> observables;
    vector<Variable> coefficients;
    offsets.push_back(constantOne);
    offsets.push_back(constantOne);

    observables.push_back(*m12);
    observables.push_back(*m13);

    coefficients.push_back(Variable("x0y0", 1.0));
    coefficients.push_back(Variable("x1y0", 0.07999, 0.01, -0.5, 0.5));
    coefficients.push_back(Variable("x2y0", -0.23732, 0.01, -0.5, 0.5));
    coefficients.push_back(Variable("x3y0", 0.10369, 0.01, -1.5, 0.5));
    coefficients.push_back(Variable("x0y1", 0.10248, 0.01, -1.5, 0.5));
    coefficients.push_back(Variable("x1y1", -0.28543, 0.01, -0.5, 0.5));
    coefficients.push_back(Variable("x2y1", 0.15058, 0.01, -1.5, 0.5));
    coefficients.push_back(Variable("x0y2", -0.20648, 0.01, -0.5, 0.5));
    coefficients.push_back(Variable("x1y2", 0.14567, 0.01, -1.5, 0.5));
    coefficients.push_back(Variable("x0y3", 0.06231, 0.01, -0.5, 0.5));

    PolynomialPdf *poly = new PolynomialPdf("efficiency", observables, coefficients, offsets, 3);
    poly->setParameterConstantness(true);

    Variable decXmin("decXmin", 6.22596);
    Variable decYmin("decYmin", 6.30722);
    Variable decZmin("decZmin", 10.82390);
    Variable conXmin("conXmin", 0.65621);
    Variable conYmin("conYmin", 0.69527);
    Variable conZmin("conZmin", 0.31764);
    Variable decXmax("decXmax", 0.79181);
    Variable decYmax("decYmax", 5.91211);
    Variable decZmax("decZmax", 1.52031);
    Variable conXmax("conXmax", 0.80918);
    Variable conYmax("conYmax", 0.22158);
    Variable conZmax("conZmax", 0.41866);

    TrigThresholdPdf *loX = new TrigThresholdPdf("loX", *m12, minDalitzX, decXmin, conXmin, false);
    TrigThresholdPdf *hiX = new TrigThresholdPdf("hiX", *m12, maxDalitzX, decXmax, conXmax, true);

    TrigThresholdPdf *loY = new TrigThresholdPdf("loY", *m13, minDalitzY, decYmin, conYmin, false);
    TrigThresholdPdf *hiY = new TrigThresholdPdf("hiY", *m13, maxDalitzY, decYmax, conYmax, true);

    TrigThresholdPdf *loZ = new TrigThresholdPdf("loZ", *m12, *m13, minDalitzZ, decZmin, conZmin, massSum, false);
    TrigThresholdPdf *hiZ = new TrigThresholdPdf("hiZ", *m12, *m13, maxDalitzZ, decZmax, conZmax, massSum, true);

    comps.clear();
    comps.push_back(poly);
    comps.push_back(loX);
    comps.push_back(hiX);
    comps.push_back(loY);
    comps.push_back(hiY);
    comps.push_back(loZ);
    comps.push_back(hiZ);
    comps.push_back(kzero_veto);
    ProdPdf *ret = new ProdPdf("efficiency_total", comps);

    // return poly;
    return ret;
}

DecayInfo3t dtop0pp{Variable("tau", 0.4101, 0.001, 0.300, 0.500),
                    Variable("xmixing", 0.0016, 0.001, 0, 0),
                    Variable("ymixing", 0.0055, 0.001, 0, 0)};

TddpPdf *makeSignalPdf(MixingTimeResolution *resolution = 0, GooPdf *eff = 0) {
    dtop0pp.motherMass   = _mD0;
    dtop0pp.daug1Mass    = piZeroMass;
    dtop0pp.daug2Mass    = piPlusMass;
    dtop0pp.daug3Mass    = piPlusMass;
    dtop0pp.meson_radius = mesonRad;

    //  dtop0pp._tau = new Variable("tau", 0.4116, 0.001, 0.300, 0.500);
    // Setting limits causes trouble with MNHESSE - why is this?
    // dtop0pp._xmixing = new Variable("xmixing", 0.01, 0.001, -0.05, 0.05);
    // dtop0pp._ymixing = new Variable("ymixing", 0.01, 0.001, -0.05, 0.05);

    // dtop0pp._tau.fixed = true;
    // dtop0pp._xmixing.fixed = true;
    // dtop0pp._ymixing.fixed = true;

    ResonancePdf *rhop = new Resonances::RBW(
        "rhop", Variable("rhop_amp_real", 1), Variable("rhop_amp_imag", 0), fixedRhoMass, fixedRhoWidth, 1, PAIR_12);

    bool fixAmps = false;

    ResonancePdf *rhom = new Resonances::RBW(
        "rhom",
        fixAmps ? Variable("rhom_amp_real", 0.714) : Variable("rhom_amp_real", 0.714, 0.001, 0, 0),
        fixAmps ? Variable("rhom_amp_imag", -0.025) : Variable("rhom_amp_imag", -0.025, 0.1, 0, 0),
        fixedRhoMass,
        fixedRhoWidth,
        1,
        PAIR_13);

    ResonancePdf *rho0
        = new Resonances::GS("rho0",
                             fixAmps ? Variable("rho0_amp_real", 0.565) : Variable("rho0_amp_real", 0.565, 0.001, 0, 0),
                             fixAmps ? Variable("rho0_amp_imag", 0.164) : Variable("rho0_amp_imag", 0.164, 0.1, 0, 0),
                             fixedRhoMass,
                             fixedRhoWidth,
                             1,
                             PAIR_23);

    Variable rho1450Mass("rhop_1450_mass", 1.465, 0.01, 1.0, 2.0);
    Variable rho1450Width("rhop_1450_width", 0.400, 0.01, 0.01, 5.0);

    ResonancePdf *rhop_1450 = new Resonances::RBW(
        "rhop_1450",
        fixAmps ? Variable("rhop_1450_amp_real", -0.174) : Variable("rhop_1450_amp_real", -0.174, 0.001, 0, 0),
        fixAmps ? Variable("rhop_1450_amp_imag", -0.117) : Variable("rhop_1450_amp_imag", -0.117, 0.1, 0, 0),
        rho1450Mass,
        rho1450Width,
        1,
        PAIR_12);

    ResonancePdf *rho0_1450 = new Resonances::RBW(
        "rho0_1450",
        fixAmps ? Variable("rho0_1450_amp_real", 0.325) : Variable("rho0_1450_amp_real", 0.325, 0.001, 0, 0),
        fixAmps ? Variable("rho0_1450_amp_imag", 0.057) : Variable("rho0_1450_amp_imag", 0.057, 0.1, 0, 0),
        rho1450Mass,
        rho1450Width,
        1,
        PAIR_23);

    ResonancePdf *rhom_1450 = new Resonances::RBW(
        "rhom_1450",
        fixAmps ? Variable("rhom_1450_amp_real", 0.788) : Variable("rhom_1450_amp_real", 0.788, 0.001, 0, 0),
        fixAmps ? Variable("rhom_1450_amp_imag", 0.226) : Variable("rhom_1450_amp_imag", 0.226, 0.1, 0, 0),
        rho1450Mass,
        rho1450Width,
        1,
        PAIR_13);

    Variable rho1700Mass("rhop_1700_mass", 1.720, 0.01, 1.6, 1.9);
    Variable rho1700Width("rhop_1700_width", 0.250, 0.01, 0.1, 1.0);

    ResonancePdf *rhop_1700 = new Resonances::RBW(
        "rhop_1700",
        fixAmps ? Variable("rhop_1700_amp_real", 2.151) : Variable("rhop_1700_amp_real", 2.151, 0.001, 0, 0),
        fixAmps ? Variable("rhop_1700_amp_imag", -0.658) : Variable("rhop_1700_amp_imag", -0.658, 0.1, 0, 0),
        rho1700Mass,
        rho1700Width,
        1,
        PAIR_12);

    ResonancePdf *rho0_1700 = new Resonances::RBW(
        "rho0_1700",
        fixAmps ? Variable("rho0_1700_amp_real", 2.400) : Variable("rho0_1700_amp_real", 2.400, 0.001, 0, 0),
        fixAmps ? Variable("rho0_1700_amp_imag", -0.734) : Variable("rho0_1700_amp_imag", -0.734, 0.1, 0, 0),
        rho1700Mass,
        rho1700Width,
        1,
        PAIR_23);

    ResonancePdf *rhom_1700 = new Resonances::RBW(
        "rhom_1700",
        fixAmps ? Variable("rhom_1700_amp_real", 1.286) : Variable("rhom_1700_amp_real", 1.286, 0.001, 0, 0),
        fixAmps ? Variable("rhom_1700_amp_imag", -1.532) : Variable("rhom_1700_amp_imag", -1.532, 0.1, 0, 0),
        rho1700Mass,
        rho1700Width,
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

    if(!drop_rho_1450) {
        dtop0pp.resonances.push_back(rhop_1450);
        dtop0pp.resonances.push_back(rho0_1450);
        dtop0pp.resonances.push_back(rhom_1450);
    }

    if(!drop_rho_1700) {
        dtop0pp.resonances.push_back(rhop_1700);
        dtop0pp.resonances.push_back(rho0_1700);
        dtop0pp.resonances.push_back(rhom_1700);
    }

    if(!drop_f0_980)
        dtop0pp.resonances.push_back(f0_980);

    if(!drop_f0_1370)
        dtop0pp.resonances.push_back(f0_1370);

    if(!drop_f0_1500)
        dtop0pp.resonances.push_back(f0_1500);

    if(!drop_f0_1710)
        dtop0pp.resonances.push_back(f0_1710);

    if(!drop_f2_1270)
        dtop0pp.resonances.push_back(f2_1270);

    if(!drop_f0_600)
        dtop0pp.resonances.push_back(f0_600);

    if(!fitMasses) {
        for(vector<ResonancePdf *>::iterator res = dtop0pp.resonances.begin(); res != dtop0pp.resonances.end(); ++res) {
            (*res)->setParameterConstantness(true);
        }
    }

    if(!eff) {
        vector<Variable> offsets;
        vector<Observable> observables;
        vector<Variable> coefficients;

        observables.push_back(*m12);
        observables.push_back(*m13);
        offsets.push_back(constantZero);
        offsets.push_back(constantZero);
        coefficients.push_back(constantOne);
        eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0);
    }

    std::vector<MixingTimeResolution *> resList;

    if(!resolution) {
        if(massd0) {
            for(int i = 0; i < mdslices; ++i) {
                sprintf(strbuffer, "coreFrac_%i", i);
                Variable coreFrac(strbuffer, 0.90, 0.001, 0.55, 0.999);
                sprintf(strbuffer, "coreBias_%i", i);
                //	Variable* coreBias = new Variable(strbuffer, 0.1, 0.001, -0.20, 0.30);
                Variable coreBias(strbuffer, -0.1, 0.001, -0.20, 0.30);
                sprintf(strbuffer, "coreScaleFactor_%i", i);
                Variable coreScaleFactor(strbuffer, 0.96, 0.001, 0.20, 1.50);
                sprintf(strbuffer, "tailScaleFactor_%i", i);
                Variable tailScaleFactor(strbuffer, 1.63, 0.001, 0.90, 6.00);
                resolution = new ThreeGaussResolution(coreFrac,
                                                      constantOne,
                                                      coreBias,
                                                      coreScaleFactor,
                                                      constantZero,
                                                      tailScaleFactor,
                                                      constantZero,
                                                      constantOne);
                resList.push_back(resolution);
            }
        } else {
            Variable coreFrac("coreFrac", 0.90, 0.001, 0.35, 0.999);
            // Variable* tailFrac = new Variable("tailFrac", 0.90, 0.001, 0.80, 1.00);
            Variable coreBias("coreBias", 0.0, 0.001, -0.20, 0.30);
            Variable coreScaleFactor("coreScaleFactor", 0.96, 0.001, 0.20, 1.50);
            // Variable* tailBias = new Variable("tailBias", 0);
            Variable tailScaleFactor("tailScaleFactor", 1.63, 0.001, 0.90, 6.00);
            // Variable* outlBias = new Variable("outlBias", 0);
            // Variable* outlScaleFactor = new Variable("outlScaleFactor", 5.0, 0.01, 0.1, 10.0);
            // resolution = new ThreeGaussResolution(coreFrac, tailFrac, coreBias, coreScaleFactor, tailBias,
            // tailScaleFactor, outlBias, outlScaleFactor);

            //      resolution = new ThreeGaussResolution(coreFrac, constantOne, coreBias, coreScaleFactor,
            //      constantZero, tailScaleFactor, constantZero, constantOne);
            //      resolution = new ThreeGaussResolution(coreFrac, constantOne, constantZero, coreScaleFactor,
            //      constantZero, tailScaleFactor, constantZero, constantOne);
            if(!doToyStudy)
                resolution = new ThreeGaussResolution(coreFrac,
                                                      constantOne,
                                                      coreBias,
                                                      coreScaleFactor,
                                                      coreBias,
                                                      tailScaleFactor,
                                                      constantZero,
                                                      constantOne);
            else {
                coreBias.setValue(0);
                coreScaleFactor.setValue(1);
                coreScaleFactor.setFixed(false);
                resolution = new ThreeGaussResolution(constantOne,
                                                      constantOne,
                                                      coreBias,
                                                      coreScaleFactor,
                                                      constantZero,
                                                      constantOne,
                                                      constantZero,
                                                      constantOne);
            }
        }
    }

    TddpPdf *mixPdf = 0;

    if(massd0)
        mixPdf = new TddpPdf("mixPdf", *dtime, *sigma, *m12, *m13, *eventNumber, dtop0pp, resList, eff, *massd0, wBkg1);
    else
        mixPdf = new TddpPdf("mixPdf", *dtime, *sigma, *m12, *m13, *eventNumber, dtop0pp, resolution, eff, wBkg1);

    return mixPdf;
}

GooPdf *makeFlatBkgDalitzPdf(bool fixem = true) {
    vector<Variable> offsets;
    vector<Observable> observables;
    vector<Variable> coefficients;
    offsets.push_back(constantZero);
    offsets.push_back(constantZero);
    observables.push_back(*m12);
    observables.push_back(*m13);
    coefficients.push_back(constantOne);

    PolynomialPdf *poly = new PolynomialPdf("flatbkgPdf", observables, coefficients, offsets, 0);
    Variable g_mean("g_mean", toyBkgTimeMean, 0.01, -0.2, 0.5);
    Variable g_sigma("g_sigma", toyBkgTimeRMS, 0.01, 0.15, 1.5);
    GooPdf *gt = new GaussianPdf("flatbkg_timepdf", *dtime, g_mean, g_sigma);
    comps.clear();
    comps.push_back(poly);
    comps.push_back(gt);
    GooPdf *ret = new ProdPdf("flatbkg_total", comps);

    if(fixem)
        ret->setParameterConstantness(true);

    return ret;
}

int runToyFit(int ifile, int nfile, bool noPlots = true) {
    if(!nfile || ifile < 0 || ifile >= 100)
        return 7; // No File or file error

    doToyStudy = true;
    //  dtime = new Variable("dtime", lowerTime, upperTime);
    dtime = new Observable("dtime", -3, 5);
    dtime->setNumBins(floor((upperTime - lowerTime) / 0.05 + 0.5));
    // dtime->getNumBins() = 200;
    // sigma = new Variable("sigma", 0, 0.8);
    sigma = new Observable("sigma", 0.099, 0.101);
    sigma->setNumBins(1);
    // Cheating way to avoid Punzi effect for toy MC. The normalisation integral is now a delta function!
    m12 = new Observable("m12", 0, 3);
    m12->setNumBins(240);
    m13 = new Observable("m13", 0, 3);
    m13->setNumBins(240);
    eventNumber = new EventNumber("eventNumber", 0, INT_MAX);
    wSig0       = new Observable("wSig0", 0, 1);

    for(int i = 0; i < nfile; i++) {
        //      sprintf(strbuffer, "dataFiles/toyPipipi0/dalitz_toyMC_%03d.txt", (i+ifile)%100);
        sprintf(strbuffer, "dataFiles/toyPipipi0/dalitz_toyMC_%03d.txt", ifile);
        toyFileName = app_ptr->get_filename(strbuffer, "examples/pipipi0DPFit");
        getToyData(toySigFraction);
    }

    // TruthResolution* dat = new TruthResolution();
    // TddpPdf* mixPdf = makeSignalPdf(dat);
    signalDalitz = makeSignalPdf();
    signalDalitz->setDataSize(data->getNumEvents(), 6); // Default 5 is fine for toys
    sig0_jsugg = new ExpPdf("sig0_jsugg", *sigma, constantZero);
    //  sig0_jsugg = makeBkg_sigma_strips(0);
    sig0_jsugg->addSpecialMask(PdfBase::ForceSeparateNorm);
    sig0_jsugg->setParameterConstantness(true);
    comps.clear();
    comps.push_back(signalDalitz);
    //  comps.push_back(sig0_jsugg);
    std::cout << "Creating overall PDF\n";
    // ProdPdf* overallSignal = new ProdPdf("overallSignal", comps);
    GooPdf *bkgPdf = makeFlatBkgDalitzPdf();
    bkgPdf->setParameterConstantness(true);

    std::vector<Observable> evtWeights;
    evtWeights.push_back(*wSig0);
    //  evtWeights.push_back(wBkg2);
    std::vector<PdfBase *> components;
    components.push_back(signalDalitz);
    components.push_back(bkgPdf);
    EventWeightedAddPdf *mixPdf = new EventWeightedAddPdf("total", evtWeights, components);
    //  GooPdf* mixPdf = overallSignal;

    mixPdf->setData(data);

    int retval;
    if(minuit1) {
        GooFit::FitManagerMinuit1 datapdf(mixPdf);
        datapdf.setMaxCalls(64000);
        datapdf.fit();
        retval = datapdf;
    } else {
        GooFit::FitManagerMinuit2 datapdf(mixPdf);
        datapdf.setMaxCalls(64000);
        datapdf.fit();
        retval = datapdf;
    }

    fmt::print("Fit results Toy fit:\n"
               "tau    : ({:.3}) fs\n"
               "xmixing: ({:.3})\%\n"
               "ymixing: ({:.3})\%\n",
               1000 * Uncertain(dtop0pp._tau),
               100 * Uncertain(dtop0pp._xmixing),
               100 * Uncertain(dtop0pp._ymixing));

    if(!noPlots)
        makeToyDalitzPlots(mixPdf);

    //  makeToyDalitzPlots(signalDalitz);
    return retval;
}

void loadDataFile(std::string fname, UnbinnedDataSet **setToFill, int effSkip) {
    if(!setToFill)
        setToFill = &data;

    std::vector<Observable> vars;
    vars.push_back(*m12);
    vars.push_back(*m13);
    vars.push_back(*dtime);
    vars.push_back(*sigma);
    vars.push_back(*eventNumber);
    vars.push_back(*wSig0);
    vars.push_back(*wBkg1);
    vars.push_back(*wBkg2);
    vars.push_back(*wBkg3);
    vars.push_back(*wBkg4);

    if(massd0)
        vars.push_back(*massd0);

    (*setToFill) = new UnbinnedDataSet(vars);
    std::ifstream reader;
    readWrapper(reader, fname.c_str());
    std::string buffer;

    while(!reader.eof()) {
        reader >> buffer;

        if(buffer == "====")
            break;

        std::cout << buffer;
    }

    double integralWeights[5] = {0, 0, 0, 0, 0};

    double dummy = 0;
    double mass  = 0;
    double delm  = 0;
    int events   = 0;

    while(!reader.eof()) {
        reader >> dummy;

        if(reader.eof())
            break;

        reader >> dummy; // m23, m(pi+ pi-), called m12 in processToyRoot convention.
        reader >> *m12;  // Already swapped according to D* charge
        reader >> *m13;

        // Errors on Dalitz variables
        reader >> dummy;
        reader >> dummy;
        reader >> dummy;

        reader >> *dtime;
        reader >> *sigma;

        if(massd0) {
            reader >> *massd0; // Md0
        } else
            reader >> mass;

        if(deltam) {
            reader >> *deltam;
        } else
            reader >> delm;

        reader >> dummy; // ProbSig
        reader >> dummy; // Dst charge
        reader >> dummy; // Run
        reader >> dummy; // Event
        reader >> *wSig0;
        reader >> *wBkg1;
        reader >> *wBkg2;
        reader >> *wBkg3;
        reader >> *wBkg4;

        if(massd0) {
            if(massd0->getValue() <= massd0->getLowerLimit())
                continue;

            if(massd0->getValue() >= massd0->getUpperLimit())
                continue;
        } else {
            // Enforce signal box on all data sets!
            if(mass <= 1.8654 + 0.0075 * md0_lower_window + md0offset)
                continue;

            if(mass >= 1.8654 + 0.0075 * md0_upper_window + md0offset)
                continue;
        }

        if(deltam) {
            if(deltam->getValue() >= deltam->getUpperLimit())
                continue;

            if(deltam->getValue() <= deltam->getLowerLimit())
                continue;
        } else {
            if(delm >= 0.1454 + 0.0003 * deltam_upper_window)
                continue;

            if(delm <= 0.1454 + 0.0003 * deltam_lower_window)
                continue;
        }

        if(dtime->getValue() < dtime->getLowerLimit())
            continue;

        if(dtime->getValue() > dtime->getUpperLimit())
            continue;

        if(sigma->getValue() > sigma->getUpperLimit())
            continue; // Lower limit is 0, and it can't be lower than that, so whatever.

        integralWeights[0] += wSig0->getValue();
        integralWeights[1] += wBkg1->getValue();
        integralWeights[2] += wBkg2->getValue();
        integralWeights[3] += wBkg3->getValue();
        integralWeights[4] += wBkg4->getValue();
        eventNumber->setValue((*setToFill)->getNumEvents());

        // See comments in TddpPdf.hh for explanation of this.
        double mistag = wSig0->getValue() + wBkg1->getValue() * luckyFrac;
        wSig0->setValue(wSig0->getValue() + wBkg1->getValue());
        wBkg1->setValue(mistag / wSig0->getValue());

        if((binEffData) && (0 == events % effSkip)) {
            double weight = weightHistogram->GetBinContent(weightHistogram->FindBin(m12->getValue(), m13->getValue()));
            // weight = 1;
            binEffData->addWeightedEvent(weight);

            // binEffData->addEvent();
            if(underlyingBins)
                underlyingBins->Fill(m12->getValue(), m13->getValue(), weight);
        } else
            (*setToFill)->addEvent();

        events++;

        if(loM23Sigma) {
            double currM23 = cpuGetM23(m12->getValue(), m13->getValue());

            if(currM23 < 1.5)
                loM23Sigma->Fill(sigma->getValue());
            else
                hiM23Sigma->Fill(sigma->getValue());

            // std::cout << "Filled sigma with " << sigma->getValue() << " " << currM23 << std::endl;
        }

        if(wSig0->getValue() > 1.0)
            std::cout << "Problem with event " << (*setToFill)->getNumEvents() << ", too-large signal weight "
                      << wSig0->getValue() << std::endl;

        /*
        if ((*setToFill)->getNumEvents() < 10) {
          std::cout << (*setToFill)->getNumEvents() << " : "
            << m12->getValue() << " "
            << m13->getValue() << " "
            << dtime->getValue() << " "
            << sigma->getValue() << " "
            << wSig0->getValue() << " "
            << wBkg1->getValue() << " "
            << wBkg2->getValue() << " "
            << wBkg3->getValue() << " "
            << wBkg4->getValue() << " "
            << std::endl;
        }
        */
    }

    std::cout << "Loaded " << (*setToFill)->getNumEvents() << " events.\n";
    std::cout << "Integrals: " << integralWeights[0] << " " << integralWeights[1] << " " << integralWeights[2] << " "
              << integralWeights[3] << " " << integralWeights[4] << "\n";
}

void makeFullFitVariables() {
    static bool exists = false;

    if(exists)
        return;

    exists = true;

    dtime = new Observable("dtime", lowerTime, upperTime);
    dtime->setNumBins(floor((upperTime - lowerTime) / 0.05 + 0.5));
    sigma = new Observable("sigma", 0, maxSigma);
    sigma->setNumBins(floor((maxSigma / 0.01) + 0.5));
    m12 = new Observable("m12", 0, 3);
    m13 = new Observable("m13", 0, 3);
    m12->setNumBins(normBinning);
    m13->setNumBins(normBinning);
    eventNumber = new EventNumber("eventNumber", 0, INT_MAX);
    wSig0       = new Observable("wSig0", 0, 1);
    wBkg1       = new Observable("wBkg1", 0, 1);
    wBkg2       = new Observable("wBkg2", 0, 1);
    wBkg3       = new Observable("wBkg3", 0, 1);
    wBkg4       = new Observable("wBkg4", 0, 1);
}

GooPdf *makeSignalJSU_gg(int idx, bool fixem = true) {
    // Values from TSigma, 'Mikhail default'.
    static int jsugg_num = -1;
    jsugg_num++;

    sprintf(strbuffer, "js_meana_%i", jsugg_num);
    Variable js_meana(strbuffer, 0.0593, 0.01, 0, 0.2);
    js_meana.setFixed(fixem);
    sprintf(strbuffer, "js_sigma_%i", jsugg_num);
    Variable js_sigma(strbuffer, 0.000474, 0.0001, 0, 0.001);
    js_sigma.setFixed(fixem);
    sprintf(strbuffer, "js_gamma_%i", jsugg_num);
    Variable js_gamma(strbuffer, -10.1942, 1, -30, 0);
    js_gamma.setFixed(fixem);
    sprintf(strbuffer, "js_delta_%i", jsugg_num);
    Variable js_delta(strbuffer, 1.4907, 0.1, 0.5, 5);
    js_delta.setFixed(fixem);
    sprintf(strbuffer, "frac_jsu_%i", jsugg_num);
    Variable frac_jsu(strbuffer, 0.9516, 0.01, 0.5, 1.0);
    frac_jsu.setFixed(fixem);
    sprintf(strbuffer, "frac_ga1_%i", jsugg_num);
    Variable frac_ga1(strbuffer, 0.001845, 0.0005, 0.0001, 0.3);
    frac_ga1.setFixed(fixem);
    sprintf(strbuffer, "g1_meana_%i", jsugg_num);
    Variable g1_meana(strbuffer, 0.2578, 0.003, 0.1, 0.5);
    g1_meana.setFixed(fixem);
    sprintf(strbuffer, "g1_sigma_%i", jsugg_num);
    Variable g1_sigma(strbuffer, 0.03086, 0.01, 0.005, 0.25);
    g1_sigma.setFixed(fixem);
    sprintf(strbuffer, "g2_meana_%i", jsugg_num);
    Variable g2_meana(strbuffer, 0.32, 0.01, 0.1, 0.5);
    g2_meana.setFixed(fixem);
    sprintf(strbuffer, "g2_sigma_%i", jsugg_num);
    Variable g2_sigma(strbuffer, 0.05825, 0.01, 0.01, 0.1);
    g2_sigma.setFixed(fixem);
    // Variable* g2_sigma = new Variable("g2_sigma", 0.5825);

    sprintf(strbuffer, "js_%i", jsugg_num);
    JohnsonSUPdf *js = new JohnsonSUPdf(strbuffer, *sigma, js_meana, js_sigma, js_gamma, js_delta);
    sprintf(strbuffer, "g1_%i", jsugg_num);
    GaussianPdf *g1 = new GaussianPdf(strbuffer, *sigma, g1_meana, g1_sigma);
    sprintf(strbuffer, "g2_%i", jsugg_num);
    // GaussianPdf*  g2 = new GaussianPdf(strbuffer, sigma, g2_meana, g2_sigma);

    weights.clear();
    weights.push_back(frac_jsu);
    // weights.push_back(frac_ga1);
    comps.clear();
    comps.push_back(js);
    comps.push_back(g1);
    // comps.push_back(g2);

    // Deal with special indices to get different starting points
    switch(idx) {
    case 1:
        g2_sigma.setUpperLimit(0.15);
        break;

    case 2:
        frac_jsu.setValue(0.80);
        break;

    // case 5:
    // g1_sigma->getLowerLimit() = 0.005;
    // break;
    case 7:
        frac_jsu.setValue(0.80);
        // frac_ga1->getValue() = 0.05;
        break;

    case 11:
        frac_ga1.setUpperLimit(0.4);
        frac_jsu.setValue(0.80);
        break;

    default:
        break;
    }

    sprintf(strbuffer, "signal_sigma_%i", jsugg_num);
    AddPdf *signal_sigma = new AddPdf(strbuffer, weights, comps);
    return signal_sigma;
}

GooPdf *makeMikhailJSU_gg(bool fixem = true) {
    // Values from TSigma, 'Mikhail default'.

    Variable *js_meana = new Variable("js_meana", 0.0593279, 0.01, 0, 0.2);
    js_meana->setFixed(fixem);
    Variable *js_sigma = new Variable("js_sigma", 0.000474171, 0.0001, 0, 0.001);
    js_sigma->setFixed(fixem);
    Variable *js_gamma = new Variable("js_gamma", -10.1942, 1, -30, 0);
    js_gamma->setFixed(fixem);
    Variable *js_delta = new Variable("js_delta", 1.4907, 0.1, 0.5, 5);
    js_delta->setFixed(fixem);
    Variable *frac_jsu = new Variable("frac_jsu", 0.951638, 0.01, 0.5, 1.0);
    frac_jsu->setFixed(fixem);
    Variable *frac_ga1 = new Variable("frac_ga1", 0.0184522, 0.00001, 0.001, 0.3);
    frac_ga1->setFixed(fixem);
    Variable *g1_meana = new Variable("g1_meana", 0.257796, 0.003, 0.1, 0.5);
    g1_meana->setFixed(fixem);
    Variable *g1_sigma = new Variable("g1_sigma", 0.0308619, 0.01, 0.005, 0.25);
    g1_sigma->setFixed(fixem);
    Variable *g2_meana = new Variable("g2_meana", 0.319952, 0.01, 0.1, 0.5);
    g2_meana->setFixed(fixem);
    Variable *g2_sigma = new Variable("g2_sigma", 0.0582524, 0.01, 0.01, 0.1);
    g2_sigma->setFixed(fixem);

    // Variable* g2_sigma = new Variable("g2_sigma", 0.5825);

    frac_ga1->setValue(frac_ga1->getValue() * (1 - frac_jsu->getValue()));

    JohnsonSUPdf *js = new JohnsonSUPdf("js", *sigma, *js_meana, *js_sigma, *js_gamma, *js_delta);
    GaussianPdf *g1  = new GaussianPdf("g1", *sigma, *g1_meana, *g1_sigma);
    GaussianPdf *g2  = new GaussianPdf("g2", *sigma, *g2_meana, *g2_sigma);

    weights.clear();
    weights.push_back(*frac_jsu);
    weights.push_back(*frac_ga1);
    comps.clear();
    comps.push_back(js);
    comps.push_back(g1);
    comps.push_back(g2);

    GooPdf *ret = new AddPdf("signal_sigma", weights, comps);
    return ret;
}

const int numSigmaBins       = 36;
TH1F **sigma_dat_hists       = 0;
TH1F **sigma_pdf_hists       = 0;
UnbinnedDataSet **sigma_data = 0;
vector<GooPdf *> jsuList;

GooPdf *makeSigmaMap() {
    sigma_dat_hists = new TH1F *[numSigmaBins];
    sigma_pdf_hists = new TH1F *[numSigmaBins];
    sigma_data      = new UnbinnedDataSet *[numSigmaBins];

    std::vector<Observable> vars;
    vars.push_back(*sigma);

    for(int i = 0; i < numSigmaBins; ++i) {
        sprintf(strbuffer, "sigma_data_hist_%i", i);
        sigma_dat_hists[i]
            = new TH1F(strbuffer, "", sigma->getNumBins(), sigma->getLowerLimit(), sigma->getUpperLimit());
        sigma_dat_hists[i]->SetStats(false);
        sigma_dat_hists[i]->SetMarkerStyle(8);
        sigma_dat_hists[i]->SetMarkerSize(0.7);
        sprintf(strbuffer, "sigma_pdf_hist_%i", i);
        sigma_pdf_hists[i]
            = new TH1F(strbuffer, "", sigma->getNumBins(), sigma->getLowerLimit(), sigma->getUpperLimit());
        sigma_pdf_hists[i]->SetStats(false);
        sigma_pdf_hists[i]->SetLineWidth(3);
        sigma_pdf_hists[i]->SetLineColor(kBlue);

        sigma_data[i] = new UnbinnedDataSet(vars);
    }

    int numEvents = data->getNumEvents();

    for(int i = 0; i < numEvents; ++i) {
        m12->setValue(data->getValue(*m12, i));
        m13->setValue(data->getValue(*m13, i));
        sigma->setValue(data->getValue(*sigma, i));

        int xbin       = (int)floor(m12->getValue() / 0.5);
        int ybin       = (int)floor(m13->getValue() / 0.5);
        int overallbin = ybin * 6 + xbin;
        sigma_dat_hists[overallbin]->Fill(sigma->getValue());
        sigma_data[overallbin]->addEvent();
    }

    // vector<GooPdf*> jsuList;
    for(int i = 0; i < numSigmaBins; ++i) {
        GooPdf *js = makeSignalJSU_gg(i, false);
        jsuList.push_back(js);

        // int xbin = i % 6;
        // int ybin = i / 6;

        /*int inDalitzPlot = 0;
        if (cpuDalitz(0.5*(xbin+0), 0.5*(ybin+0), _mD0, piZeroMass, piPlusMass, piPlusMass)) inDalitzPlot++;
        if (cpuDalitz(0.5*(xbin+1), 0.5*(ybin+0), _mD0, piZeroMass, piPlusMass, piPlusMass)) inDalitzPlot++;
        if (cpuDalitz(0.5*(xbin+0), 0.5*(ybin+1), _mD0, piZeroMass, piPlusMass, piPlusMass)) inDalitzPlot++;
        if (cpuDalitz(0.5*(xbin+1), 0.5*(ybin+1), _mD0, piZeroMass, piPlusMass, piPlusMass)) inDalitzPlot++;
        */
        if(0 == sigma_data[i]->getNumEvents())
            js->setParameterConstantness(true);
        else {
            std::cout << "\n\nAbout to start fit of sigma box " << i << std::endl;
            js->setData(sigma_data[i]);
            if(minuit1) {
                GooFit::FitManagerMinuit1 currpdf(js);
                currpdf.fit();
            } else {
                GooFit::FitManagerMinuit2 currpdf(js);
                currpdf.fit();
            }
            js->setParameterConstantness(true);
            // js->clearCurrentFit();
            std::cout << "Done with sigma box " << i << "\n";
        }
    }

    vector<Observable> obses;
    obses.push_back(*m12);
    obses.push_back(*m13);
    vector<double> limits;
    limits.push_back(0);
    limits.push_back(0);
    vector<double> binSizes;
    binSizes.push_back(0.5);
    binSizes.push_back(0.5);
    vector<int> numBins;
    numBins.push_back(6);
    numBins.push_back(6);
    BinTransformPdf *mapFunction = new BinTransformPdf("mapFunction", obses, limits, binSizes, numBins);

    return new MappedPdf("sigmaMap", mapFunction, jsuList);
}

GooPdf *make1BinSigmaMap() {
    sigma_dat_hists = new TH1F *[1];
    sigma_pdf_hists = new TH1F *[1];
    sigma_data      = new UnbinnedDataSet *[1];

    std::vector<Observable> vars;
    vars.push_back(*sigma);

    for(int i = 0; i < 1; ++i) {
        sprintf(strbuffer, "sigma_data_hist_%i", i);
        sigma_dat_hists[i]
            = new TH1F(strbuffer, "", sigma->getNumBins(), sigma->getLowerLimit(), sigma->getUpperLimit());
        sigma_dat_hists[i]->SetStats(false);
        sigma_dat_hists[i]->SetMarkerStyle(8);
        sigma_dat_hists[i]->SetMarkerSize(0.7);
        sprintf(strbuffer, "sigma_pdf_hist_%i", i);
        sigma_pdf_hists[i]
            = new TH1F(strbuffer, "", sigma->getNumBins(), sigma->getLowerLimit(), sigma->getUpperLimit());
        sigma_pdf_hists[i]->SetStats(false);
        sigma_pdf_hists[i]->SetLineWidth(3);
        sigma_pdf_hists[i]->SetLineColor(kBlue);

        sigma_data[i] = new UnbinnedDataSet(vars);
    }

    int numEvents = data->getNumEvents();

    for(int i = 0; i < numEvents; ++i) {
        m12->setValue(data->getValue(*m12, i));
        m13->setValue(data->getValue(*m13, i));
        sigma->setValue(data->getValue(*sigma, i));

        int overallbin = 0;
        sigma_dat_hists[overallbin]->Fill(sigma->getValue());
        sigma_data[overallbin]->addEvent();
    }

    // vector<GooPdf*> jsuList;
    for(int i = 0; i < 1; ++i) {
        GooPdf *js = makeSignalJSU_gg(i, false);
        jsuList.push_back(js);

        if(0 == sigma_data[i]->getNumEvents())
            js->setParameterConstantness(true);
        else {
            std::cout << "\n\nAbout to start fit of sigma box " << i << std::endl;
            js->setData(sigma_data[i]);
            if(minuit1) {
                GooFit::FitManagerMinuit1 currpdf(js);
                currpdf.fit();
            } else {
                GooFit::FitManagerMinuit2 currpdf(js);
                currpdf.fit();
            }
            js->setParameterConstantness(true);
            // js->clearCurrentFit();
            std::cout << "Done with sigma box " << i << "\n";
        }
    }

    vector<Observable> obses;
    obses.push_back(*m12);
    obses.push_back(*m13);
    vector<double> limits;
    limits.push_back(0);
    limits.push_back(0);
    vector<double> binSizes;
    binSizes.push_back(3);
    binSizes.push_back(3);
    vector<int> numBins;
    numBins.push_back(1);
    numBins.push_back(1);
    BinTransformPdf *mapFunction = new BinTransformPdf("mapFunction", obses, limits, binSizes, numBins);

    return new MappedPdf("sigmaMap", mapFunction, jsuList);
}

GooPdf *make4BinSigmaMap() {
    sigma_dat_hists = new TH1F *[4];
    sigma_pdf_hists = new TH1F *[4];
    sigma_data      = new UnbinnedDataSet *[4];

    std::vector<Observable> vars;
    vars.push_back(*sigma);

    for(int i = 0; i < 4; ++i) {
        sprintf(strbuffer, "sigma_data_hist_%i", i);
        sigma_dat_hists[i]
            = new TH1F(strbuffer, "", sigma->getNumBins(), sigma->getLowerLimit(), sigma->getUpperLimit());
        sigma_dat_hists[i]->SetStats(false);
        sigma_dat_hists[i]->SetMarkerStyle(8);
        sigma_dat_hists[i]->SetMarkerSize(0.7);
        sprintf(strbuffer, "sigma_pdf_hist_%i", i);
        sigma_pdf_hists[i]
            = new TH1F(strbuffer, "", sigma->getNumBins(), sigma->getLowerLimit(), sigma->getUpperLimit());
        sigma_pdf_hists[i]->SetStats(false);
        sigma_pdf_hists[i]->SetLineWidth(3);
        sigma_pdf_hists[i]->SetLineColor(kBlue);

        sigma_data[i] = new UnbinnedDataSet(vars);
    }

    int numEvents = data->getNumEvents();

    for(int i = 0; i < numEvents; ++i) {
        m12->setValue(data->getValue(*m12, i));
        m13->setValue(data->getValue(*m13, i));
        sigma->setValue(data->getValue(*sigma, i));

        int xbin       = (int)floor(m12->getValue() / 1.5);
        int ybin       = (int)floor(m13->getValue() / 1.5);
        int overallbin = ybin * 2 + xbin;
        sigma_dat_hists[overallbin]->Fill(sigma->getValue());
        sigma_data[overallbin]->addEvent();
    }

    // vector<GooPdf*> jsuList;
    for(int i = 0; i < 4; ++i) {
        GooPdf *js = makeSignalJSU_gg(i, false);
        jsuList.push_back(js);

        if(0 == sigma_data[i]->getNumEvents())
            js->setParameterConstantness(true);
        else {
            std::cout << "\n\nAbout to start fit of sigma box " << i << std::endl;
            js->setData(sigma_data[i]);
            if(minuit1) {
                GooFit::FitManagerMinuit1 currpdf(js);
                currpdf.fit();
            } else {
                GooFit::FitManagerMinuit2 currpdf(js);
                currpdf.fit();
            }
            js->setParameterConstantness(true);
            // js->clearCurrentFit();
            std::cout << "Done with sigma box " << i << "\n";
        }
    }

    vector<Observable> obses;
    obses.push_back(*m12);
    obses.push_back(*m13);
    vector<double> limits;
    limits.push_back(0);
    limits.push_back(0);
    vector<double> binSizes;
    binSizes.push_back(1.5);
    binSizes.push_back(1.5);
    vector<int> numBins;
    numBins.push_back(2);
    numBins.push_back(2);
    BinTransformPdf *mapFunction = new BinTransformPdf("mapFunction", obses, limits, binSizes, numBins);

    return new MappedPdf("sigmaMap", mapFunction, jsuList);
}

void coarseBin(TH2F &dalitzHist, int grain) {
    // Move from too-fine-to-see binning down to reasonable binning in Dalitz plots.
    for(int i = 1; i < m12->getNumBins(); i += grain) {
        for(int j = 1; j < m13->getNumBins(); j += grain) {
            double total = 0;

            for(int k = 0; k < grain; ++k) {
                for(int l = 0; l < grain; ++l) {
                    total += dalitzHist.GetBinContent(i + k, j + l);
                }
            }

            for(int k = 0; k < grain; ++k) {
                for(int l = 0; l < grain; ++l) {
                    dalitzHist.SetBinContent(i + k, j + l, total);
                }
            }
        }
    }
}

struct BigBin {
    int xbin;
    int ybin;
    int width;
    int height;
    double getContent(TH2F *plot);
};

double BigBin::getContent(TH2F *plot) {
    double ret = 0;

    // std::cout << "getContent with " << width << " " << height << " " << xbin << " " << ybin <<std::endl;
    for(unsigned int i = 0; i < width; ++i) {
        for(unsigned int j = 0; j < height; ++j) {
            // std::cout << i << ", " << j << std::endl;
            if(xbin + i > plot->GetNbinsX())
                continue;

            if(ybin + j > plot->GetNbinsY())
                continue;

            ret += plot->GetBinContent(xbin + i, ybin + j);
        }
    }

    // std::cout << "Total " << ret << std::endl;
    return ret;
}

struct ChisqInfo {
    ChisqInfo();
    double chisq;
    int dof;
    TH2F *contribPlot;
};

ChisqInfo::ChisqInfo()
    : chisq(0)
    , dof(0)
    , contribPlot(0) {}

ChisqInfo *getAdaptiveChisquare(TH2F *datPlot, TH2F *pdfPlot) {
    bool acceptable = false;
    int binSize     = 1;
    vector<BigBin> binlist;
    double limit = 26;

    while(!acceptable) {
        binlist.clear();
        std::cout << "Attempting bin generation with size " << binSize << std::endl;

        for(int xbin = 1; xbin <= datPlot->GetNbinsX(); xbin += binSize) {
            for(int ybin = 1; ybin <= datPlot->GetNbinsY(); ybin += binSize) {
                double lox  = datPlot->GetXaxis()->GetBinLowEdge(xbin + 0);
                double hix  = datPlot->GetXaxis()->GetBinLowEdge(xbin + 1 + binSize);
                double loy  = datPlot->GetYaxis()->GetBinLowEdge(ybin + 0);
                double hiy  = datPlot->GetYaxis()->GetBinLowEdge(ybin + 1 + binSize);
                bool corner = false;

                if(cpuDalitz(lox, loy, _mD0, piZeroMass, piPlusMass, piPlusMass))
                    corner = true;
                else if(cpuDalitz(lox, hiy, _mD0, piZeroMass, piPlusMass, piPlusMass))
                    corner = true;
                else if(cpuDalitz(hix, loy, _mD0, piZeroMass, piPlusMass, piPlusMass))
                    corner = true;
                else if(cpuDalitz(hix, hiy, _mD0, piZeroMass, piPlusMass, piPlusMass))
                    corner = true;

                if(!corner)
                    continue;

                BigBin curr;
                curr.xbin  = xbin;
                curr.ybin  = ybin;
                curr.width = curr.height = binSize;
                binlist.push_back(curr);
            }
        }

        acceptable = true;

        for(vector<BigBin>::iterator bin = binlist.begin(); bin != binlist.end(); ++bin) {
            if((*bin).getContent(datPlot) >= limit)
                continue;

            acceptable = false;
            binSize *= 2;
            std::cout << "Couldn't get good bins, retry.\n";
            break;
        }
    }

    std::cout << "Good bins at size " << binSize << ", beginning splits.\n";

    // Now attempt to split bins.
    int numSplits = 1;

    while(0 < numSplits) {
        numSplits = 0;
        vector<BigBin> newbins;

        for(vector<BigBin>::iterator bin = binlist.begin(); bin != binlist.end(); ++bin) {
            if(1 == (*bin).width * (*bin).height) {
                newbins.push_back(*bin);
                continue;
            }

            BigBin lolef;
            BigBin lorig;
            BigBin hilef;
            BigBin hirig;
            lolef.xbin = (*bin).xbin;
            hilef.xbin = (*bin).xbin;
            lorig.xbin = (*bin).xbin + (*bin).width / 2;
            hirig.xbin = (*bin).xbin + (*bin).width / 2;
            lolef.ybin = (*bin).ybin;
            hilef.ybin = (*bin).ybin + (*bin).height / 2;
            lorig.ybin = (*bin).ybin;
            hirig.ybin = (*bin).ybin + (*bin).height / 2;

            lolef.width  = (*bin).width / 2;
            lorig.width  = (*bin).width / 2;
            hilef.width  = (*bin).width / 2;
            hirig.width  = (*bin).width / 2;
            lolef.height = (*bin).height / 2;
            lorig.height = (*bin).height / 2;
            hilef.height = (*bin).height / 2;
            hirig.height = (*bin).height / 2;

            int mask = 0;

            if(limit < lolef.getContent(datPlot))
                mask += 1;

            if(limit < lorig.getContent(datPlot))
                mask += 2;

            if(limit < hilef.getContent(datPlot))
                mask += 4;

            if(limit < hirig.getContent(datPlot))
                mask += 8;

            if(mask != 15) {
                newbins.push_back(*bin);
            } else {
                newbins.push_back(lolef);
                newbins.push_back(lorig);
                newbins.push_back(hilef);
                newbins.push_back(hirig);
                numSplits++;
            }
        }

        binlist.clear();

        for(vector<BigBin>::iterator i = newbins.begin(); i != newbins.end(); ++i)
            binlist.push_back(*i);

        std::cout << "Split " << numSplits << " bins.\n";
    }

    ChisqInfo *ret   = new ChisqInfo();
    ret->dof         = binlist.size();
    ret->contribPlot = new TH2F("contribPlot",
                                "",
                                datPlot->GetNbinsX(),
                                datPlot->GetXaxis()->GetBinLowEdge(1),
                                datPlot->GetXaxis()->GetBinLowEdge(datPlot->GetNbinsX() + 1),
                                datPlot->GetNbinsY(),
                                datPlot->GetYaxis()->GetBinLowEdge(1),
                                datPlot->GetYaxis()->GetBinLowEdge(datPlot->GetNbinsY() + 1));

    double totalDat = 0;
    double totalPdf = 0;

    for(vector<BigBin>::iterator bin = binlist.begin(); bin != binlist.end(); ++bin) {
        double dat  = (*bin).getContent(datPlot);
        double pdf  = (*bin).getContent(pdfPlot);
        double term = (dat - pdf) / sqrt(dat);
        ret->chisq += term * term;

        /*
        std::cout << "Bin (" << (*bin).xbin << ", " << (*bin).ybin << ") "
              << (*bin).width << " " << (*bin).height << " : "
              << dat << " " << pdf << " "
              << term << std::endl;
        */
        for(int i = 0; i < (*bin).width; ++i) {
            for(int j = 0; j < (*bin).height; ++j) {
                bool corner = false;
                double lox  = datPlot->GetXaxis()->GetBinLowEdge((*bin).xbin + i);
                double hix  = datPlot->GetXaxis()->GetBinLowEdge((*bin).xbin + i + 1);
                double loy  = datPlot->GetYaxis()->GetBinLowEdge((*bin).ybin + j);
                double hiy  = datPlot->GetYaxis()->GetBinLowEdge((*bin).ybin + j + 1);

                if(cpuDalitz(lox, loy, _mD0, piZeroMass, piPlusMass, piPlusMass))
                    corner = true;
                else if(cpuDalitz(lox, hiy, _mD0, piZeroMass, piPlusMass, piPlusMass))
                    corner = true;
                else if(cpuDalitz(hix, loy, _mD0, piZeroMass, piPlusMass, piPlusMass))
                    corner = true;
                else if(cpuDalitz(hix, hiy, _mD0, piZeroMass, piPlusMass, piPlusMass))
                    corner = true;

                if(!corner)
                    continue;

                ret->contribPlot->SetBinContent((*bin).xbin + i, (*bin).ybin + j, term);
            }
        }

        totalPdf += pdf;
        totalDat += dat;
    }

    return ret;
}

void makeToyDalitzPlots(GooPdf *overallSignal, std::string plotdir) {
    std::string call = "mkdir -p " + plotdir;
    system(call.c_str());

    foo->cd();

    TH1F dtime_dat_hist("dtime_dat_hist", "", dtime->getNumBins(), dtime->getLowerLimit(), dtime->getUpperLimit());
    dtime_dat_hist.SetStats(false);
    dtime_dat_hist.SetMarkerStyle(8);
    dtime_dat_hist.SetMarkerSize(1.2);
    dtime_dat_hist.GetXaxis()->SetTitle("Decay time [ps]");
    dtime_dat_hist.GetYaxis()->SetTitle("Events / 50 fs");
    TH1F dtime_pdf_hist("dtime_pdf_hist", "", dtime->getNumBins(), dtime->getLowerLimit(), dtime->getUpperLimit());
    dtime_pdf_hist.SetStats(false);
    dtime_pdf_hist.SetLineColor(kBlue);
    dtime_pdf_hist.SetLineWidth(3);
    TH1F dtime_sig_hist("dtime_sig_hist", "", dtime->getNumBins(), dtime->getLowerLimit(), dtime->getUpperLimit());
    dtime_sig_hist.SetStats(false);
    dtime_sig_hist.SetLineColor(kRed);
    dtime_sig_hist.SetLineWidth(3);
    TH1F dtime_bg_hist("dtime_bg_hist", "", dtime->getNumBins(), dtime->getLowerLimit(), dtime->getUpperLimit());
    dtime_bg_hist.SetStats(false);
    dtime_bg_hist.SetLineColor(kMagenta);
    dtime_bg_hist.SetLineWidth(3);

    double totalPdf     = 0;
    double totalPdf_sig = 0;
    double totalPdf_bg  = 0;
    double totalDat     = 0;
    double totalSigProb = 0;
    double totalBGProb  = 0;

    for(unsigned int evt = 0; evt < data->getNumEvents(); ++evt) {
        double currTime = data->getValue(*dtime, evt);
        dtime_dat_hist.Fill(currTime);
        totalSigProb += data->getValue(*wSig0, evt);
        totalBGProb += 1 - data->getValue(*wSig0, evt);
        totalDat++;
    }

    std::cout << "totalData = " << totalDat << ", totalSigProb = " << totalSigProb << std::endl;
    std::vector<Observable> vars;
    vars.push_back(*m12);
    vars.push_back(*m13);
    vars.push_back(*dtime);
    vars.push_back(*sigma);
    vars.push_back(*eventNumber);
    vars.push_back(*wSig0);
    UnbinnedDataSet currData(vars);
    sigma->setValue(0.1);
    wSig0->setValue(totalSigProb / totalDat);
    int evtCounter = 0;

    for(int i = 0; i < m12->getNumBins(); ++i) {
        m12->setValue(m12->getLowerLimit()
                      + (m12->getUpperLimit() - m12->getLowerLimit()) * (i + 0.5) / m12->getNumBins());

        for(int j = 0; j < m13->getNumBins(); ++j) {
            m13->setValue(m13->getLowerLimit()
                          + (m13->getUpperLimit() - m13->getLowerLimit()) * (j + 0.5) / m13->getNumBins());

            if(!cpuDalitz(m12->getValue(), m13->getValue(), _mD0, piZeroMass, piPlusMass, piPlusMass))
                continue;

            for(int l = 0; l < dtime->getNumBins(); ++l) {
                dtime->setValue(dtime->getLowerLimit()
                                + (dtime->getUpperLimit() - dtime->getLowerLimit()) * (l + 0.5) / dtime->getNumBins());
                eventNumber->setValue(evtCounter);
                evtCounter++;
                currData.addEvent();
            }
        }
    }

    GOOFIT_INFO("Adding {} signal events from toy", currData.getNumEvents());

    overallSignal->setData(&currData);
    signalDalitz->setDataSize(currData.getNumEvents(), 6);
    std::vector<std::vector<double>> pdfValues = overallSignal->getCompProbsAtDataPoints();

    for(unsigned int j = 0; j < pdfValues[0].size(); ++j) {
        double currTime = currData.getValue(*dtime, j);
        dtime_sig_hist.Fill(currTime, pdfValues[1][j]);
        dtime_bg_hist.Fill(currTime, pdfValues[2][j]);
        totalPdf += pdfValues[0][j];
        totalPdf_sig += pdfValues[1][j];
        totalPdf_bg += pdfValues[2][j];
    }

    for(int i = 1; i <= dtime->getNumBins(); ++i) {
        dtime_sig_hist.SetBinContent(i, dtime_sig_hist.GetBinContent(i) * totalSigProb / totalPdf_sig);
        dtime_bg_hist.SetBinContent(i, dtime_bg_hist.GetBinContent(i) * totalBGProb / totalPdf_bg);
        dtime_pdf_hist.SetBinContent(i, dtime_sig_hist.GetBinContent(i) + dtime_bg_hist.GetBinContent(i));
    }

    foo->cd();
    dtime_dat_hist.Draw("p");
    dtime_pdf_hist.Draw("lsame");
    dtime_bg_hist.SetLineStyle(2);
    dtime_bg_hist.Draw("lsame");
    dtime_sig_hist.SetLineStyle(3);
    dtime_sig_hist.Draw("lsame");

    foo->SaveAs((plotdir + "/dtime_fit.png").c_str());
    foo->SetLogy(true);
    foo->SaveAs((plotdir + "/dtime_fit_log.png").c_str());
    foo->SetLogy(false);
}

void makeDalitzPlots(GooPdf *overallSignal, std::string plotdir = "./plots_from_mixfit/") {
    std::string mkplotdir{"mkdir " + plotdir};
    system(mkplotdir.c_str());
    foo->cd();

    TH1F dtime_dat_hist("dtime_dat_hist", "", dtime->getNumBins(), dtime->getLowerLimit(), dtime->getUpperLimit());
    dtime_dat_hist.SetStats(false);
    dtime_dat_hist.SetMarkerStyle(8);
    dtime_dat_hist.SetMarkerSize(1.2);
    dtime_dat_hist.GetXaxis()->SetTitle("Decay time [ps]");
    dtime_dat_hist.GetYaxis()->SetTitle("Events / 50 fs");
    TH1F dtime_pdf_hist("dtime_pdf_hist", "", dtime->getNumBins(), dtime->getLowerLimit(), dtime->getUpperLimit());
    dtime_pdf_hist.SetStats(false);
    dtime_pdf_hist.SetLineColor(kBlue);
    dtime_pdf_hist.SetLineWidth(3);

    TH1F sigma_dat_hist("sigma_dat_hist", "", sigma->getNumBins(), sigma->getLowerLimit(), sigma->getUpperLimit());
    sigma_dat_hist.SetStats(false);
    sigma_dat_hist.SetMarkerStyle(8);
    sigma_dat_hist.SetMarkerSize(1.2);
    sigma_dat_hist.GetXaxis()->SetTitle("Decay time error [ps]");
    sigma_dat_hist.GetYaxis()->SetTitle("Events / 8 fs");
    TH1F sigma_pdf_hist("sigma_pdf_hist", "", sigma->getNumBins(), sigma->getLowerLimit(), sigma->getUpperLimit());
    sigma_pdf_hist.SetStats(false);
    sigma_pdf_hist.SetLineColor(kBlue);
    sigma_pdf_hist.SetLineWidth(3);

    TH1F m12_dat_hist("m12_dat_hist", "", m12->getNumBins(), m12->getLowerLimit(), m12->getUpperLimit());
    m12_dat_hist.SetStats(false);
    m12_dat_hist.SetMarkerStyle(8);
    m12_dat_hist.SetMarkerSize(1.2);
    m12_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{0}) [GeV]");
    m12_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    TH1F m12_pdf_hist("m12_pdf_hist", "", m12->getNumBins(), m12->getLowerLimit(), m12->getUpperLimit());
    m12_pdf_hist.SetStats(false);
    m12_pdf_hist.SetLineColor(kBlue);
    m12_pdf_hist.SetLineWidth(3);

    TH1F m13_dat_hist("m13_dat_hist", "", m13->getNumBins(), m13->getLowerLimit(), m13->getUpperLimit());
    m13_dat_hist.SetStats(false);
    m13_dat_hist.SetMarkerStyle(8);
    m13_dat_hist.SetMarkerSize(1.2);
    m13_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{0}) [GeV]");
    m13_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    TH1F m13_pdf_hist("m13_pdf_hist", "", m13->getNumBins(), m13->getLowerLimit(), m13->getUpperLimit());
    m13_pdf_hist.SetStats(false);
    m13_pdf_hist.SetLineColor(kBlue);
    m13_pdf_hist.SetLineWidth(3);

    TH1F m23_dat_hist("m23_dat_hist", "", m13->getNumBins(), m13->getLowerLimit(), m13->getUpperLimit());
    m23_dat_hist.SetStats(false);
    m23_dat_hist.SetMarkerStyle(8);
    m23_dat_hist.SetMarkerSize(1.2);
    m23_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{+}) [GeV]");
    m23_dat_hist.GetYaxis()->SetTitle("Events / 12.5 MeV");
    TH1F m23_pdf_hist("m23_pdf_hist", "", m13->getNumBins(), m13->getLowerLimit(), m13->getUpperLimit());
    m23_pdf_hist.SetStats(false);
    m23_pdf_hist.SetLineColor(kBlue);
    m23_pdf_hist.SetLineWidth(3);

    TH2F dalitzpm_dat_hist("dalitzpm_dat_hist",
                           "",
                           m12->getNumBins(),
                           m12->getLowerLimit(),
                           m12->getUpperLimit(),
                           m13->getNumBins(),
                           m13->getLowerLimit(),
                           m13->getUpperLimit());
    dalitzpm_dat_hist.SetStats(false);
    dalitzpm_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{0}) [GeV]");
    dalitzpm_dat_hist.GetYaxis()->SetTitle("m^{2}(#pi^{-} #pi^{0}) [GeV]");
    TH2F dalitzpm_pdf_hist("dalitzpm_pdf_hist",
                           "",
                           m12->getNumBins(),
                           m12->getLowerLimit(),
                           m12->getUpperLimit(),
                           m13->getNumBins(),
                           m13->getLowerLimit(),
                           m13->getUpperLimit());
    dalitzpm_pdf_hist.SetStats(false);

    TH2F dalitzp0_dat_hist("dalitzp0_dat_hist",
                           "",
                           m12->getNumBins(),
                           m12->getLowerLimit(),
                           m12->getUpperLimit(),
                           m13->getNumBins(),
                           m13->getLowerLimit(),
                           m13->getUpperLimit());
    dalitzp0_dat_hist.SetStats(false);
    dalitzp0_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{0}) [GeV]");
    dalitzp0_dat_hist.GetYaxis()->SetTitle("m^{2}(#pi^{-} #pi^{+}) [GeV]");
    TH2F dalitzp0_pdf_hist("dalitzp0_pdf_hist",
                           "",
                           m12->getNumBins(),
                           m12->getLowerLimit(),
                           m12->getUpperLimit(),
                           m13->getNumBins(),
                           m13->getLowerLimit(),
                           m13->getUpperLimit());
    dalitzp0_pdf_hist.SetStats(false);

    TH2F dalitzm0_dat_hist("dalitzm0_dat_hist",
                           "",
                           m12->getNumBins(),
                           m12->getLowerLimit(),
                           m12->getUpperLimit(),
                           m13->getNumBins(),
                           m13->getLowerLimit(),
                           m13->getUpperLimit());
    dalitzm0_dat_hist.SetStats(false);
    dalitzm0_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{0}) [GeV]");
    dalitzm0_dat_hist.GetYaxis()->SetTitle("m^{2}(#pi^{+} #pi^{-}) [GeV]");
    TH2F dalitzm0_pdf_hist("dalitzm0_pdf_hist",
                           "",
                           m12->getNumBins(),
                           m12->getLowerLimit(),
                           m12->getUpperLimit(),
                           m13->getNumBins(),
                           m13->getLowerLimit(),
                           m13->getUpperLimit());
    dalitzm0_pdf_hist.SetStats(false);

    TH1F *bkg3_pdfs[6];
    TH1F *bkg3_data[6];
    double num_sigma_dat[6];
    double num_sigma_pdf[6];

    for(int i = 0; i < 6; ++i) {
        sprintf(strbuffer, "bkg_sigma_pdf_%i", i);
        bkg3_pdfs[i] = new TH1F(strbuffer, "", sigma->getNumBins(), sigma->getLowerLimit(), sigma->getUpperLimit());
        sprintf(strbuffer, "bkg_sigma_dat_%i", i);
        bkg3_data[i] = new TH1F(strbuffer, "", sigma->getNumBins(), sigma->getLowerLimit(), sigma->getUpperLimit());

        num_sigma_dat[i] = 0;
        num_sigma_pdf[i] = 0;

        bkg3_data[i]->SetStats(false);
        bkg3_data[i]->SetMarkerStyle(8);
        bkg3_data[i]->SetMarkerSize(1.2);
        bkg3_data[i]->GetXaxis()->SetTitle("Decay time error [ps]");
        bkg3_data[i]->GetYaxis()->SetTitle("Events / 8 fs");

        bkg3_pdfs[i]->SetStats(false);
        bkg3_pdfs[i]->SetLineColor(kBlue);
        bkg3_pdfs[i]->SetLineWidth(3);
    }

    double totalPdf = 0;
    double totalDat = 0;

    for(unsigned int evt = 0; evt < data->getNumEvents(); ++evt) {
        double currTime = data->getValue(*dtime, evt);
        dtime_dat_hist.Fill(currTime);

        double currSigma = data->getValue(*sigma, evt);
        sigma_dat_hist.Fill(currSigma);

        double currm12 = data->getValue(*m12, evt);
        m12_dat_hist.Fill(currm12);

        double currm13 = data->getValue(*m13, evt);
        m13_dat_hist.Fill(currm13);

        dalitzpm_dat_hist.Fill(currm12, currm13);

        double currm23 = cpuGetM23(currm12, currm13);
        m23_dat_hist.Fill(currm23);

        dalitzp0_dat_hist.Fill(currm12, currm23);
        dalitzm0_dat_hist.Fill(currm13, currm23);

        totalDat++;

        int m23bin = (int)floor(currm23 / 0.5);
        bkg3_data[m23bin]->Fill(currSigma);
        num_sigma_dat[m23bin]++;
    }

    double maxBinContent = 0;
    int bestI            = 0;
    int bestJ            = 0;

    for(int i = 1; i <= m12->getNumBins(); ++i) {
        for(int j = 1; j <= m13->getNumBins(); ++j) {
            double curr = dalitzpm_dat_hist.GetBinContent(i, j);

            if(curr < maxBinContent)
                continue;

            maxBinContent = curr;
            bestI         = i;
            bestJ         = j;
        }
    }

    std::cout << "Max bin content: " << maxBinContent << " (" << bestI << ", " << bestJ << ")\n";

    bool dependsOnSigma           = true;
    std::vector<Observable> obses = overallSignal->getObservables();

    if(std::find(obses.begin(), obses.end(), *sigma) == obses.end())
        dependsOnSigma = false;

    // overallSignal->setDebugMask(1);

    wBkg1->setValue(0);
    const int division = 2;

    for(int half = 0; half < division; ++half) {
        std::vector<Observable> vars;
        vars.push_back(*m12);
        vars.push_back(*m13);
        vars.push_back(*dtime);
        vars.push_back(*sigma);
        vars.push_back(*eventNumber);
        vars.push_back(*wBkg1);

        UnbinnedDataSet currData(vars);

        sigma->setValue(0.5);
        int evtCounter = 0;

        for(int i = 0; i < m12->getNumBins(); ++i) {
            m12->setValue(m12->getLowerLimit()
                          + (m12->getUpperLimit() - m12->getLowerLimit()) * (i + 0.5) / m12->getNumBins());

            for(int j = 0; j < m13->getNumBins(); ++j) {
                m13->setValue(m13->getLowerLimit()
                              + (m13->getUpperLimit() - m13->getLowerLimit()) * (j + 0.5) / m13->getNumBins());

                if(!cpuDalitz(m12->getValue(), m13->getValue(), _mD0, piZeroMass, piPlusMass, piPlusMass))
                    continue;

                for(int l = half; l < dtime->getNumBins(); l += division) {
                    dtime->setValue(dtime->getLowerLimit()
                                    + (dtime->getUpperLimit() - dtime->getLowerLimit()) * (l + 0.5)
                                          / dtime->getNumBins());
                    eventNumber->setValue(evtCounter);
                    evtCounter++;
                    currData.addEvent();
                }
            }
        }

        for(int k = 0; k < sigma->getNumBins(); ++k) {
            sigma->setValue(sigma->getLowerLimit()
                            + (sigma->getUpperLimit() - sigma->getLowerLimit()) * (k + 0.5) / sigma->getNumBins());

            if(0 == k % 10)
                std::cout << "sigma iteration " << half << " " << k << std::endl;

            currData.setValueForAllEvents(*sigma);
            overallSignal->setData(&currData);

            if(0 == k) {
                if(signalDalitz)
                    signalDalitz->setDataSize(currData.getNumEvents(), 6);

                if(incsum1)
                    incsum1->setDataSize(currData.getNumEvents(), 6);

                if(incsum2)
                    incsum2->setDataSize(currData.getNumEvents(), 6);

                if(incsum3)
                    incsum3->setDataSize(currData.getNumEvents(), 6);

                if(incsum4)
                    incsum4->setDataSize(currData.getNumEvents(), 6);

                if(incsum5)
                    incsum5->setDataSize(currData.getNumEvents(), 6);

                if(incsum6)
                    incsum6->setDataSize(currData.getNumEvents(), 6);
            }

            std::vector<std::vector<double>> pdfValues = overallSignal->getCompProbsAtDataPoints();

            for(unsigned int j = 0; j < pdfValues[0].size(); ++j) {
                double currTime = currData.getValue(*dtime, j);
                dtime_pdf_hist.Fill(currTime, pdfValues[0][j]);

                double currSigma = currData.getValue(*sigma, j);
                sigma_pdf_hist.Fill(currSigma, pdfValues[0][j]);

                // Um... these two are switched? Weirdness...
                double currm12 = currData.getValue(*m13, j);
                m12_pdf_hist.Fill(currm12, pdfValues[0][j]);

                double currm13 = currData.getValue(*m12, j);
                m13_pdf_hist.Fill(currm13, pdfValues[0][j]);
                dalitzpm_pdf_hist.Fill(currm12, currm13, pdfValues[0][j]);

                double currm23 = cpuGetM23(currm12, currm13);
                m23_pdf_hist.Fill(currm23, pdfValues[0][j]);
                dalitzp0_pdf_hist.Fill(currm12, currm23, pdfValues[0][j]);
                dalitzm0_pdf_hist.Fill(currm13, currm23, pdfValues[0][j]);

                int currM23Bin = (int)(currm23 / 0.5);
                bkg3_pdfs[currM23Bin]->Fill(currSigma, pdfValues[0][j]);
                num_sigma_pdf[currM23Bin] += pdfValues[0][j];

                totalPdf += pdfValues[0][j];

                if(std::isnan(pdfValues[0][j])) {
                    std::cout << "Major problem: " << k << " " << j << std::endl;
                    assert(false);
                }

                if(std::isinf(pdfValues[0][j])) {
                    std::cout << "Infinity " << k << " " << j << std::endl;
                    assert(false);
                }
            }

            // If PDF doesn't depend on sigma, don't project from that dimension.
            if(!dependsOnSigma)
                break;
        }
    }

    // std::cout << "Final values: " << totalDat << " " << totalPdf << std::endl;

    for(int i = 1; i <= dtime->getNumBins(); ++i) {
        dtime_pdf_hist.SetBinContent(i, dtime_pdf_hist.GetBinContent(i) * totalDat / totalPdf);
    }

    for(int i = 1; i <= sigma->getNumBins(); ++i) {
        sigma_pdf_hist.SetBinContent(i, sigma_pdf_hist.GetBinContent(i) * totalDat / totalPdf);

        for(int j = 0; j < 6; ++j) {
            bkg3_pdfs[j]->SetBinContent(i, bkg3_pdfs[j]->GetBinContent(i) * num_sigma_dat[j] / num_sigma_pdf[j]);
        }
    }

    for(int i = 1; i <= m12->getNumBins(); ++i) {
        m12_pdf_hist.SetBinContent(i, m12_pdf_hist.GetBinContent(i) * totalDat / totalPdf);
    }

    for(int i = 1; i <= m13->getNumBins(); ++i) {
        m13_pdf_hist.SetBinContent(i, m13_pdf_hist.GetBinContent(i) * totalDat / totalPdf);
        m23_pdf_hist.SetBinContent(i, m23_pdf_hist.GetBinContent(i) * totalDat / totalPdf);
    }

    for(int i = 1; i <= m12->getNumBins(); ++i) {
        for(int j = 1; j <= m13->getNumBins(); ++j) {
            dalitzpm_pdf_hist.SetBinContent(i, j, dalitzpm_pdf_hist.GetBinContent(i, j) * totalDat / totalPdf);
            dalitzp0_pdf_hist.SetBinContent(i, j, dalitzp0_pdf_hist.GetBinContent(i, j) * totalDat / totalPdf);
            dalitzm0_pdf_hist.SetBinContent(i, j, dalitzm0_pdf_hist.GetBinContent(i, j) * totalDat / totalPdf);
        }
    }

    ChisqInfo *chisq = getAdaptiveChisquare(&dalitzpm_dat_hist, &dalitzpm_pdf_hist);

    std::cout << "Chisquare: " << chisq->chisq << " / " << chisq->dof << std::endl;
    foodal->cd();
    foodal->SetLogz(false);
    chisq->contribPlot->SetStats(false);
    chisq->contribPlot->Draw("colz");
    foodal->SaveAs((plotdir + "/chisq.png").c_str());
    foo->cd();

    coarseBin(dalitzpm_pdf_hist, 2);
    coarseBin(dalitzpm_dat_hist, 2);
    coarseBin(dalitzp0_pdf_hist, 2);
    coarseBin(dalitzp0_dat_hist, 2);
    coarseBin(dalitzm0_pdf_hist, 2);
    coarseBin(dalitzm0_dat_hist, 2);

    dtime_dat_hist.Draw("p");
    dtime_pdf_hist.Draw("lsame");
    foo->SaveAs((plotdir + "/dtime_fit.png").c_str());
    foo->SetLogy(true);
    foo->SaveAs((plotdir + "/dtime_fit_log.png").c_str());
    foo->SetLogy(false);

    for(int i = 0; i < 6; ++i) {
        if(!dependsOnSigma) {
            bkg3_data[i]->Draw("p");
        } else if(bkg3_data[i]->GetMaximum() > bkg3_pdfs[i]->GetMaximum()) {
            bkg3_data[i]->Draw("p");
            bkg3_pdfs[i]->Draw("lsame");
        } else {
            bkg3_pdfs[i]->Draw("l");
            bkg3_data[i]->Draw("psame");
        }

        sprintf(strbuffer, "%i", i);
        TText slicenum;
        slicenum.DrawTextNDC(0.2, 0.8, strbuffer);

        foo->SaveAs((plotdir + bkg3_pdfs[i]->GetName() + ".png").c_str());
        foo->SetLogy(true);
        foo->SaveAs((plotdir + bkg3_pdfs[i]->GetName() + "_log.png").c_str());
        foo->SetLogy(false);
    }

    m13_dat_hist.Draw("p");
    m13_pdf_hist.Draw("lsame");
    foo->SaveAs((plotdir + "/m13_fit.png").c_str());
    foo->SetLogy(true);
    foo->SaveAs((plotdir + "/m13_fit_log.png").c_str());
    foo->SetLogy(false);

    sigma_dat_hist.Draw("p");
    sigma_pdf_hist.Draw("lsame");
    foo->SaveAs((plotdir + "/sigma_fit.png").c_str());
    foo->SetLogy(true);
    foo->SaveAs((plotdir + "/sigma_fit_log.png").c_str());
    foo->SetLogy(false);

    m12_dat_hist.Draw("p");
    m12_pdf_hist.Draw("lsame");
    foo->SaveAs((plotdir + "/m12_fit.png").c_str());
    foo->SetLogy(true);
    foo->SaveAs((plotdir + "/m12_fit_log.png").c_str());
    foo->SetLogy(false);

    m23_dat_hist.Draw("p");
    m23_pdf_hist.Draw("lsame");
    foo->SaveAs((plotdir + "/m23_fit.png").c_str());
    foo->SetLogy(true);
    foo->SaveAs((plotdir + "/m23_fit_log.png").c_str());
    foo->SetLogy(false);

    foodal->cd();
    dalitzpm_dat_hist.Draw("colz");

    for(int slice = 0; slice < 6; ++slice) {
        double line_m12 = cpuGetM23(0, 0.5 * (slice + 1));
        TLine sliceLine;
        sliceLine.SetLineWidth(2);
        sliceLine.DrawLine(0, line_m12, line_m12, 0);
        sprintf(strbuffer, "%i", slice);
        TText sliceNum;
        sliceNum.DrawText(0.25, line_m12 - 0.25, strbuffer);
    }

    foodal->SaveAs((plotdir + "/dalitzpm_dat.png").c_str());
    foodal->SetLogz(true);
    foodal->SaveAs((plotdir + "/dalitzpm_dat_log.png").c_str());
    foodal->SetLogz(false);
    dalitzp0_dat_hist.Draw("colz");
    foodal->SaveAs((plotdir + "/dalitzp0_dat.png").c_str());
    foodal->SetLogz(true);
    foodal->SaveAs((plotdir + "/dalitzp0_dat_log.png").c_str());
    foodal->SetLogz(false);
    dalitzm0_dat_hist.Draw("colz");
    foodal->SaveAs((plotdir + "/dalitzm0_dat.png").c_str());
    foodal->SetLogz(true);
    foodal->SaveAs((plotdir + "/dalitzm0_dat_log.png").c_str());
    foodal->SetLogz(false);

    dalitzpm_pdf_hist.Draw("colz");
    foodal->SaveAs((plotdir + "/dalitzpm_pdf.png").c_str());
    foodal->SetLogz(true);
    foodal->SaveAs((plotdir + "/dalitzpm_pdf_log.png").c_str());
    foodal->SetLogz(false);
    dalitzpm_pdf_hist.GetZaxis()->SetRangeUser(0, dalitzpm_dat_hist.GetMaximum());
    dalitzpm_pdf_hist.Draw("colz");
    foodal->SaveAs((plotdir + "/dalitzpm_pdf_matched.png").c_str());
    foodal->SetLogz(true);
    foodal->SaveAs((plotdir + "/dalitzpm_pdf_matched_log.png").c_str());
    foodal->SetLogz(false);
    dalitzp0_pdf_hist.GetZaxis()->SetRangeUser(0, dalitzp0_dat_hist.GetMaximum());
    dalitzp0_pdf_hist.Draw("colz");
    foodal->SaveAs((plotdir + "/dalitzp0_pdf.png").c_str());
    foodal->SetLogz(true);
    foodal->SaveAs((plotdir + "/dalitzp0_pdf_log.png").c_str());
    foodal->SetLogz(false);
    dalitzm0_pdf_hist.GetZaxis()->SetRangeUser(0, dalitzm0_dat_hist.GetMaximum());
    dalitzm0_pdf_hist.Draw("colz");
    foodal->SaveAs((plotdir + "/dalitzm0_pdf.png").c_str());
    foodal->SetLogz(true);
    foodal->SaveAs((plotdir + "/dalitzm0_pdf_log.png").c_str());
    foodal->SetLogz(false);

    TH1F pull_pm_hist("pull_pm_hist", "", 100, -5, 5);
    pull_pm_hist.GetXaxis()->SetTitle("(Data - PDF) / sqrt(Data)");
    pull_pm_hist.GetYaxis()->SetTitle("Bins / 0.1");

    for(int i = 1; i <= m12->getNumBins(); ++i) {
        double m12loedge
            = m12->getLowerLimit() + ((m12->getUpperLimit() - m12->getLowerLimit()) / m12->getNumBins()) * (i - 1);
        double m12hiedge
            = m12->getLowerLimit() + ((m12->getUpperLimit() - m12->getLowerLimit()) / m12->getNumBins()) * (i);

        for(int j = 1; j <= m13->getNumBins(); ++j) {
            double m13loedge
                = m13->getLowerLimit() + ((m13->getUpperLimit() - m13->getLowerLimit()) / m13->getNumBins()) * (j - 1);

            if(!cpuDalitz(m12loedge, m13loedge, _mD0, piZeroMass, piPlusMass, piPlusMass)) {
                dalitzpm_dat_hist.SetBinContent(i, j, 0);
                continue;
            }

            if(!cpuDalitz(m12hiedge, m13loedge, _mD0, piZeroMass, piPlusMass, piPlusMass)) {
                dalitzpm_dat_hist.SetBinContent(i, j, 0);
                continue;
            }

            double m13hiedge
                = m13->getLowerLimit() + ((m13->getUpperLimit() - m13->getLowerLimit()) / m13->getNumBins()) * (j);

            if(!cpuDalitz(m12loedge, m13hiedge, _mD0, piZeroMass, piPlusMass, piPlusMass)) {
                dalitzpm_dat_hist.SetBinContent(i, j, 0);
                continue;
            }

            if(!cpuDalitz(m12hiedge, m13hiedge, _mD0, piZeroMass, piPlusMass, piPlusMass)) {
                dalitzpm_dat_hist.SetBinContent(i, j, 0);
                continue;
            }

            double dat = dalitzpm_dat_hist.GetBinContent(i, j);
            double pdf = dalitzpm_pdf_hist.GetBinContent(i, j);

            double pullval = (dat - pdf) / sqrt(max(1.0, dat));
            dalitzpm_dat_hist.SetBinContent(i, j, pullval);
            pull_pm_hist.Fill(pullval);
        }

        for(int j = 1; j <= m13->getNumBins(); ++j) {
            double m23loedge
                = m13->getLowerLimit() + ((m13->getUpperLimit() - m13->getLowerLimit()) / m13->getNumBins()) * (j - 1);

            // To get 12, 23 instead of 12, 13, just exchange 1 and 2.
            if(!cpuDalitz(m12loedge, m23loedge, _mD0, piPlusMass, piZeroMass, piPlusMass)) {
                dalitzp0_dat_hist.SetBinContent(i, j, 0);
                continue;
            }

            if(!cpuDalitz(m12hiedge, m23loedge, _mD0, piPlusMass, piZeroMass, piPlusMass)) {
                dalitzp0_dat_hist.SetBinContent(i, j, 0);
                continue;
            }

            double m23hiedge
                = m13->getLowerLimit() + ((m13->getUpperLimit() - m13->getLowerLimit()) / m13->getNumBins()) * (j);

            if(!cpuDalitz(m12loedge, m23hiedge, _mD0, piPlusMass, piZeroMass, piPlusMass)) {
                dalitzp0_dat_hist.SetBinContent(i, j, 0);
                continue;
            }

            if(!cpuDalitz(m12hiedge, m23hiedge, _mD0, piPlusMass, piZeroMass, piPlusMass)) {
                dalitzp0_dat_hist.SetBinContent(i, j, 0);
                continue;
            }

            double dat = dalitzp0_dat_hist.GetBinContent(i, j);
            double pdf = dalitzp0_pdf_hist.GetBinContent(i, j);

            dalitzp0_dat_hist.SetBinContent(i, j, (dat - pdf) / sqrt(max(1.0, dat)));
        }

        // NB, this exploits symmetry 12 and 13 by treating the outer loop as 13.
        for(int j = 1; j <= m13->getNumBins(); ++j) {
            double m23loedge
                = m13->getLowerLimit() + ((m13->getUpperLimit() - m13->getLowerLimit()) / m13->getNumBins()) * (j - 1);

            if(!cpuDalitz(m12loedge, m23loedge, _mD0, piPlusMass, piZeroMass, piPlusMass)) {
                dalitzm0_dat_hist.SetBinContent(i, j, 0);
                continue;
            }

            if(!cpuDalitz(m12hiedge, m23loedge, _mD0, piPlusMass, piZeroMass, piPlusMass)) {
                dalitzm0_dat_hist.SetBinContent(i, j, 0);
                continue;
            }

            double m23hiedge
                = m13->getLowerLimit() + ((m13->getUpperLimit() - m13->getLowerLimit()) / m13->getNumBins()) * (j);

            if(!cpuDalitz(m12loedge, m23hiedge, _mD0, piPlusMass, piZeroMass, piPlusMass)) {
                dalitzm0_dat_hist.SetBinContent(i, j, 0);
                continue;
            }

            if(!cpuDalitz(m12hiedge, m23hiedge, _mD0, piPlusMass, piZeroMass, piPlusMass)) {
                dalitzm0_dat_hist.SetBinContent(i, j, 0);
                continue;
            }

            double dat = dalitzm0_dat_hist.GetBinContent(i, j);
            double pdf = dalitzm0_pdf_hist.GetBinContent(i, j);

            dalitzm0_dat_hist.SetBinContent(i, j, (dat - pdf) / sqrt(max(1.0, dat)));
        }
    }

    dalitzpm_dat_hist.GetZaxis()->SetRangeUser(-5, 5);
    dalitzpm_dat_hist.Draw("colz");
    foodal->SaveAs((plotdir + "/dalitzpm_pull.png").c_str());
    dalitzp0_dat_hist.GetZaxis()->SetRangeUser(-5, 5);
    dalitzp0_dat_hist.Draw("colz");
    foodal->SaveAs((plotdir + "/dalitzp0_pull.png").c_str());
    dalitzm0_dat_hist.GetZaxis()->SetRangeUser(-5, 5);
    dalitzm0_dat_hist.Draw("colz");
    foodal->SaveAs((plotdir + "/dalitzm0_pull.png").c_str());

    foo->cd();

    pull_pm_hist.Draw("");
    foo->SaveAs((plotdir + "/pull_pm_hist.png").c_str());
}

GooPdf *make_m23_transform() {
    // This is complicated. I want to make a function such that the parameters
    // depend on position in the Dalitz plot. I want to use stripes of m23, because
    // that seems to be the real shape of the dependence, and I haven't got so much
    // data for bkg3; so stripes work better than boxes of m12, m13. So I need a
    // transform from m12, m13 to m23. Then I need a transform from m23 to a bin number.
    // Finally I need a transform from bin number to function. Keep the tongue straight
    // in the mouth, now!

    vector<Observable> obses;
    vector<Variable> offsets;
    vector<Variable> coefficients;
    vector<PdfBase *> components;
    vector<double> limits;
    vector<double> binSizes;
    vector<int> numBins;

    // First the transform to m23.
    // m23 = _mD02 + piZeroMass*piZeroMass + piPlusMass*piPlusMass + piPlusMass*piPlusMass - massPZ - massPM

    obses.push_back(*m12);
    obses.push_back(*m13);
    coefficients.push_back(constantBigM);
    coefficients.push_back(constantMinusOne);
    coefficients.push_back(constantMinusOne);
    offsets.push_back(constantZero);
    offsets.push_back(constantZero);
    PolynomialPdf *m23_transform = new PolynomialPdf("m23_transform", obses, coefficients, offsets, 1);

    // Now create the BinTransform which will make bins out of m23 values.
    obses.clear();
    coefficients.clear();
    obses.push_back(*m12); // Fake dependence; CompositePdf will create a fake event using this index.
    limits.push_back(0);   // Bins of m23 start at 0.
    binSizes.push_back(3.0 / m23Slices);
    numBins.push_back(m23Slices);
    BinTransformPdf *m23_binMap = new BinTransformPdf("m23_binMap", obses, limits, binSizes, numBins);

    // Now make a composite, so that the BinTransform takes the Polynomial result as input.
    CompositePdf *m23_composite = new CompositePdf("m23_composite", m23_transform, m23_binMap);
    return m23_composite;
}

GooPdf *makeSigmaHists() {
    GooPdf *m23_composite = make_m23_transform();
    std::vector<std::unique_ptr<BinnedDataSet>> sigmaHists;

    for(int i = 0; i < m23Slices; ++i)
        sigmaHists.emplace_back(new BinnedDataSet(*sigma));

    std::ifstream reader;
    std::string fname = app_ptr->get_filename("./dataFiles/signalMC_truth_mm_0.txt", "examples/pipipi0DPFit");
    readWrapper(reader, fname);
    std::string buffer;

    while(!reader.eof()) {
        reader >> buffer;

        if(buffer == "====")
            break;

        std::cout << buffer;
    }

    double dummy = 0;
    double m23   = 0;

    while(!reader.eof()) {
        reader >> dummy;
        reader >> m23;
        reader >> dummy;
        reader >> dummy;
        reader >> dummy;
        reader >> dummy;
        reader >> dummy;
        reader >> *dtime;
        reader >> *sigma;

        reader >> dummy;
        reader >> dummy;
        reader >> dummy;
        reader >> dummy;
        reader >> dummy;
        reader >> dummy;
        reader >> dummy;
        reader >> dummy;
        reader >> dummy;
        reader >> dummy;
        reader >> dummy;

        if(dtime->getValue() < dtime->getLowerLimit())
            continue;

        if(dtime->getValue() > dtime->getUpperLimit())
            continue;

        int bin = (int)floor((m23 / 3.0) * m23Slices);
        sigmaHists[bin]->addEvent();
    }

    vector<GooPdf *> jsuList;

    for(int i = 0; i < m23Slices; ++i) {
        sprintf(strbuffer, "signal_sigma_hist_%i", i);
        SmoothHistogramPdf *hist = new SmoothHistogramPdf(strbuffer, sigmaHists[i].get(), constantZero);
        jsuList.push_back(hist);
    }

    return new MappedPdf("signalSigmaHist", m23_composite, jsuList);
}

GooPdf *makeBkg_sigma_strips(int bkgnum) {
    GooPdf *m23_composite = make_m23_transform();

    vector<GooPdf *> jsuList;
    vector<ConvolutionPdf *> convList;
    bool useShare = false;

    for(int i = 0; i < m23Slices; ++i) {
        sprintf(strbuffer, "bkg%i_sigma_slice%i_expalpha", bkgnum, i);
        Variable exp_alpha(strbuffer, 7.50, 0.10, 0, 10.00);
        sprintf(strbuffer, "bkg%i_sigma_slice%i_gauss_meana", bkgnum, i);
        Variable g_meana(strbuffer, 0.20, 0.01, 0.00, 0.80);
        sprintf(strbuffer, "bkg%i_sigma_slice%i_gauss_sigma", bkgnum, i);
        Variable g_sigma(strbuffer, 0.03, 0.01, 0.01, 0.20);

        sprintf(strbuffer, "bkg%i_sigma_slice%i_conv", bkgnum, i);
        ExpGausPdf *expfunc = new ExpGausPdf(strbuffer, *sigma, g_meana, g_sigma, exp_alpha);
        jsuList.push_back(expfunc);
    }

    if(useShare) {
        for(vector<ConvolutionPdf *>::iterator conv = convList.begin(); conv != convList.end(); ++conv) {
            (*conv)->registerOthers(convList);
        }
    }

    sprintf(strbuffer, "bkg%i_sigma_map", bkgnum);
    return new MappedPdf(strbuffer, m23_composite, jsuList);
}

void createWeightHistogram() {
    weightHistogram = new TH2F("weightHistogram",
                               "",
                               m12->getNumBins(),
                               m12->getLowerLimit(),
                               m12->getUpperLimit(),
                               m13->getNumBins(),
                               m13->getLowerLimit(),
                               m13->getUpperLimit());
    weightHistogram->SetStats(false);
    double step12 = (m12->getUpperLimit() - m12->getLowerLimit()) / m12->getNumBins();
    double step13 = (m13->getUpperLimit() - m13->getLowerLimit()) / m13->getNumBins();

    for(int i = 1; i <= m12->getNumBins(); ++i) {
        for(int j = 1; j < m13->getNumBins(); ++j) {
            double maxCount = 0;
            double count    = 0;

            for(double currM12 = m12->getLowerLimit() + step12 * (i - 1) + 0.05 * step12;
                currM12 < m12->getLowerLimit() + step12 * i;
                currM12 += 0.1 * step12) {
                for(double currM13 = m13->getLowerLimit() + step13 * (j - 1) + 0.05 * step13;
                    currM13 < m13->getLowerLimit() + step13 * j;
                    currM13 += 0.1 * step13) {
                    maxCount++;

                    if(!cpuDalitz(currM12, currM13, _mD0, piZeroMass, piPlusMass, piPlusMass))
                        continue;

                    count++;
                }
            }

            if(0.1 > maxCount)
                continue;

            count /= maxCount;
            weightHistogram->SetBinContent(i, j, count);
        }
    }

    underlyingBins = new TH2F("underlyingBins",
                              "",
                              m12->getNumBins(),
                              m12->getLowerLimit(),
                              m12->getUpperLimit(),
                              m13->getNumBins(),
                              m13->getLowerLimit(),
                              m13->getUpperLimit());
    underlyingBins->SetStats(false);
}

GooPdf *makeOverallSignal() {
    makeKzeroVeto();

    int oldBins1 = m12->getNumBins();
    int oldBins2 = m13->getNumBins();
    // Too fine a binning here leads to bad results due to fluctuations.
    m12->setNumBins(120);
    m13->setNumBins(120);
    vector<Observable> lvars;
    lvars.push_back(*m12);
    lvars.push_back(*m13);
    binEffData = new BinnedDataSet(lvars);
    createWeightHistogram();
    std::cout << "Loading efficiency data" << std::endl;
    std::string fname = app_ptr->get_filename("./dataFiles/efficiency_flat.txt", "examples/pipipi0DPFit");
    loadDataFile(fname, &effdata);

    if(saveEffPlot) {
        system("mkdir plots_from_mixfit");
        foodal->cd();
        underlyingBins->Draw("colz");
        foodal->SaveAs("plots_from_mixfit/efficiency_bins.png");
        foodal->SetLogz(true);
        foodal->SaveAs("plots_from_mixfit/efficiency_bins_log.png");
        foo->cd();
    }

    GooPdf *eff = 0;

    // Polynomial version
    if(polyEff)
        eff = makeEfficiencyPoly();
    // SmoothHistogram version
    else
        eff = makeEfficiencyPdf();

    eff->setData(effdata);

    if(minuit1) {
        GooFit::FitManagerMinuit1 effpdf(eff);
        effpdf.fit();
    } else {
        GooFit::FitManagerMinuit2 effpdf(eff);
        effpdf.fit();
    }

    eff->setParameterConstantness(true);
    delete binEffData;
    binEffData = nullptr;
    delete effdata;
    effdata = nullptr;

    m12->setNumBins(oldBins1);
    m13->setNumBins(oldBins2);

    comps.clear();
    comps.push_back(eff);
    comps.push_back(kzero_veto);
    ProdPdf *effWithVeto = new ProdPdf("effWithVeto", comps);

    std::cout << "Creating signal PDF\n";
    signalDalitz = makeSignalPdf(0, effWithVeto);

    std::cout << "Creating sigma PDF\n";

    // sig0_jsugg = makeSigmaMap();
    // sig0_jsugg = make1BinSigmaMap();
    // sig0_jsugg = make4BinSigmaMap();
    // sig0_jsugg = makeMikhailJSU_gg();
    if(useHistogramSigma)
        sig0_jsugg = makeSigmaHists();
    else
        sig0_jsugg = makeBkg_sigma_strips(0);

    sig0_jsugg->addSpecialMask(PdfBase::ForceSeparateNorm);
    // sig0_jsugg = makeSignalJSU_gg(-1, false);

    /*
    sig0_jsugg->setData(data);
    FitManager jsupdf(sig0_jsugg);
    gettimeofday(&startTime, nullptr);
    jsupdf.fit();
    gettimeofday(&stopTime, nullptr);
    timersub(&stopTime, &startTime, &totalTime);
    std::cout << "Time for sigma fit : " << totalTime.tv_sec + totalTime.tv_usec/1000000.0 << " seconds." << std::endl;
    */
    sprintf(strbuffer, "signal_sigma_%islices_pdf.txt", m23Slices);
    fname = app_ptr->get_filename(strbuffer, "examples/pipipi0DPFit");
    readFromFile(sig0_jsugg, strbuffer);
    sig0_jsugg->setParameterConstantness(true);

    comps.clear();
    comps.push_back(signalDalitz);
    comps.push_back(sig0_jsugg);
    std::cout << "Creating overall PDF\n";
    ProdPdf *overallSignal = new ProdPdf("overallSignal", comps);
    return overallSignal;
}

int runTruthMCFit(std::string fname, bool noPlots = true) {
    makeFullFitVariables();

    std::cout << "Loading MC data from " << fname << std::endl;
    loadDataFile(fname);
    GooPdf *overallSignal = makeOverallSignal();

    signalDalitz->setDataSize(data->getNumEvents()); // Default 5 is ok here, no event weighting
    overallSignal->setData(data);
    // overallSignal->setDebugMask(1);

    int retval;
    if(minuit1) {
        GooFit::FitManagerMinuit1 datapdf(overallSignal);
        datapdf.fit();
        retval = datapdf;
    } else {
        GooFit::FitManagerMinuit2 datapdf(overallSignal);
        datapdf.fit();
        retval = datapdf;
    }

    // overallSignal->setDebugMask(0);

    fmt::print("Fit results Toy fit TruthMC fit:\n"
               "tau    : {:.3}\n"
               "xmixing: ({:.3})\%\n"
               "ymixing: ({:.3})\%\n",
               Uncertain(dtop0pp._tau),
               100 * Uncertain(dtop0pp._xmixing),
               100 * Uncertain(dtop0pp._ymixing));

    if(noPlots)
        return retval;

    makeDalitzPlots(overallSignal, "./plots_from_mixfit/fullMCfit/");
    return retval;
}

int runGeneratedMCFit(std::string fname, int genResolutions, double dplotres) {
    makeFullFitVariables();
    std::cout << "Loading (generated) MC data from " << fname << std::endl;
    dtime->setUpperLimit(6);
    dtime->setLowerLimit(0);
    sigma->setUpperLimit(1.0);
    loadDataFile(fname);

    TRandom donram(42);
    std::vector<Observable> vars = data->getObservables();
    UnbinnedDataSet *smearedData = new UnbinnedDataSet(vars);

    if(0 != genResolutions) {
        int numEvents = data->getNumEvents();

        for(int i = 0; i < numEvents; ++i) {
            data->loadEvent(i);

            double smear1 = 0;
            double smear2 = 0;

            if(DplotRes & genResolutions) {
                if(100 > smearedData->getNumEvents())
                    std::cout << "Before smear: " << m12->getValue() << " " << m13->getValue();

                smear1 = donram.Gaus(0, dplotres);
                smear2 = donram.Gaus(0, dplotres);
            }

            // if (cpuDalitz(m12->getValue() + smear1, m13->getValue() + smear2, _mD0, piZeroMass, piPlusMass,
            // piPlusMass)) {
            m12->setValue(m12->getValue() + smear1);
            m13->setValue(m13->getValue() + smear2);

            //}
            if(100 > smearedData->getNumEvents())
                std::cout << " After smear: " << m12->getValue() << " " << m13->getValue() << "\n";

            smearedData->addEvent();
        }

        // delete data;
        // data = smearedData;
    } else
        smearedData = data;

    /*
    vector<Variable*> lvars;
    lvars.push_back(m12);
    lvars.push_back(m13);
    binEffData = new BinnedDataSet(lvars);

    TH2F genEff("genEff", "", m12->getNumBins(), m12->getLowerLimit(), m12->getUpperLimit(), m13->getNumBins(),
    m13->getLowerLimit(), m13->getUpperLimit());
    TH2F resEff("resEff", "", m12->getNumBins(), m12->getLowerLimit(), m12->getUpperLimit(), m13->getNumBins(),
    m13->getLowerLimit(), m13->getUpperLimit());

    double xstep = (m12->getUpperLimit() - m12->getLowerLimit()) / m12->getNumBins();
    double ystep = (m13->getUpperLimit() - m13->getLowerLimit()) / m13->getNumBins();
    for (int i = 0; i < m12->getNumBins(); ++i) {
      double lox = m12->getLowerLimit() + i*xstep;
      double hix = m12->getLowerLimit() + (i+1)*xstep;
      if (0 == i%10) std::cout << "Generating efficiency for " << i << std::endl;
      for (int j = 0; j < m13->getNumBins(); ++j) {
        double loy = m13->getLowerLimit() + j*ystep;
        double hiy = m13->getLowerLimit() + (j+1)*ystep;

        bool corner = false;
        if (cpuDalitz(lox, loy, _mD0, piZeroMass, piPlusMass, piPlusMass)) corner = true;
        else if (cpuDalitz(lox, hiy, _mD0, piZeroMass, piPlusMass, piPlusMass)) corner = true;
        else if (cpuDalitz(hix, loy, _mD0, piZeroMass, piPlusMass, piPlusMass)) corner = true;
        else if (cpuDalitz(hix, hiy, _mD0, piZeroMass, piPlusMass, piPlusMass)) corner = true;
        if (!corner) continue;

        for (int k = 0; k < 10000; ++k) {
    double genx = donram.Uniform(lox, hix);
    double geny = donram.Uniform(loy, hiy);
    if (!cpuDalitz(genx, geny, _mD0, piZeroMass, piPlusMass, piPlusMass)) continue;
    genEff.Fill(genx, geny);
    if (DplotRes & genResolutions) {
      genx += donram.Gaus(0, dplotres);
      geny += donram.Gaus(0, dplotres);
    }
    if (!cpuDalitz(genx, geny, _mD0, piZeroMass, piPlusMass, piPlusMass)) continue;
    resEff.Fill(genx, geny);
        }
      }
    }

    for (int i = 0; i < m12->getNumBins(); ++i) {
      for (int j = 0; j < m13->getNumBins(); ++j) {
        double gen = genEff.GetBinContent(i+1, j+1);
        if (0.1 > gen) continue;
        double res = resEff.GetBinContent(i+1, j+1);
        m12->getValue() = m12->getLowerLimit() + (i+0.5)*xstep;
        m13->getValue() = m13->getLowerLimit() + (j+0.5)*ystep;
        binEffData->setBinContent(binEffData->getBinNumber(), (res/gen));
        resEff.SetBinContent(i+1, j+1, (res/gen));
      }
    }

    Variable* effSmoothing = new Variable("effSmoothing", 0.0);
    SmoothHistogramPdf* eff = new SmoothHistogramPdf("efficiency", binEffData, effSmoothing);
    foodal->cd();
    resEff.SetStats(false);
    foodal->SetLogz(true);
    resEff.GetZaxis()->SetRangeUser(0.99, 1);
    resEff.Draw("colz");
    foodal->SaveAs("gen_efficiency.png");
    foodal->SetLogz(false);
    foo->cd();
    */

    int oldBins1 = m12->getNumBins();
    int oldBins2 = m13->getNumBins();
    m12->setNumBins(120);
    m13->setNumBins(120);

    createWeightHistogram();
    weightHistogram->Draw("colz");
    foodal->SaveAs("./plots_from_mixfit/efficiency_weights.png");

    vector<Observable> lvars;
    lvars.push_back(*m12);
    lvars.push_back(*m13);
    binEffData = new BinnedDataSet(lvars);
    fname      = app_ptr->get_filename("dataFiles/efficiency_gen.txt", "examples/pipipi0DPFit");
    loadDataFile(fname, &effdata, 1);
    GooPdf *eff = makeEfficiencyPdf();
    m12->setNumBins(oldBins1);
    m13->setNumBins(oldBins2);

    // eff->setData(effdata);
    // FitManager effpdf(eff);
    // effpdf.fit();
    // eff->setParameterConstantness(true);
    // binEffData = 0;
    // delete effdata; effdata = 0;

    TruthResolution *res = new TruthResolution();
    signalDalitz         = makeSignalPdf(res, eff);
    signalDalitz->setDataSize(smearedData->getNumEvents()); // Default 5 is ok here, no event weighting
    signalDalitz->setData(smearedData);

    /*
    std::vector<std::vector<double> > pdfValues1;
    signalDalitz->getCompProbsAtDataPoints(pdfValues1);

    signalDalitz->setDataSize(smearedData->getNumEvents()); // Default 5 is ok here, no event weighting
    signalDalitz->setData(smearedData);
    std::vector<std::vector<double> > pdfValues2;
    signalDalitz->getCompProbsAtDataPoints(pdfValues2);

    double nll1 = 0;
    double nll2 = 0;
    for (int i = 0; i < data->getNumEvents(); ++i) {
      nll1 += log(pdfValues1[0][i]);
      nll2 += log(pdfValues2[0][i]);
      data->loadEvent(i);
      if ((100 > i) || (fabs(pdfValues1[0][i] - pdfValues2[0][i]) > 0.5)) {
        double eff1 = binEffData->getBinContent(binEffData->getBinNumber());
        std::cout << i << ": " << m12->getValue() << " " << m13->getValue() << " -> " << pdfValues1[0][i] << " " << eff1
    << " " << binEffData->getBinNumber();
        smearedData->loadEvent(i);
        eff1 = binEffData->getBinContent(binEffData->getBinNumber());
        std::cout <<     " | " << m12->getValue() << " " << m13->getValue() << " -> " << pdfValues2[0][i] << " " << eff1
    << " " << binEffData->getBinNumber();
        std::cout << std::endl;
      }
    }
    std::cout << "Final NLLs: " << nll1  << " " << nll2 << std::endl;
    */

    int retval;
    if(minuit1) {
        GooFit::FitManagerMinuit1 datapdf(signalDalitz);
        datapdf.setMaxCalls(64000);
        datapdf.fit();
        retval = datapdf;
    } else {
        GooFit::FitManagerMinuit2 datapdf(signalDalitz);
        datapdf.setMaxCalls(64000);
        datapdf.fit();
        retval = datapdf;
    }

    fmt::print("Fit results Canonical fit:\n"
               "tau    : {:.3}\n"
               "xmixing: ({:.3})\%\n"
               "ymixing: ({:.3})\%\n",
               Uncertain(dtop0pp._tau),
               100 * Uncertain(dtop0pp._xmixing),
               100 * Uncertain(dtop0pp._ymixing));

    // All this relies on exact formatting of the input data files; it's fragile.
    double inputx              = 1;
    double inputy              = 1;
    std::string::size_type pos = fname.find("mm");

    if(pos != std::string::npos)
        inputx = inputy = -1;
    else {
        pos = fname.find("mp");

        if(pos != std::string::npos)
            inputx = -1;
        else {
            pos = fname.find("pm");

            if(pos != std::string::npos)
                inputy = -1;
            else {
                pos = fname.find("pp");
                assert(pos != std::string::npos);
            }
        }
    }

    std::string ident = fname.substr(pos, 4);
    sprintf(strbuffer, "result_%s_%f", ident.c_str(), dplotres);
    ofstream writer;
    writer.open(strbuffer);
    writer << inputx << " " << 100 * dtop0pp._xmixing.getValue() << " " << 100 * dtop0pp._xmixing.getError() << " "
           << inputy << " " << 100 * dtop0pp._ymixing.getValue() << " " << 100 * dtop0pp._ymixing.getError()
           << std::endl;
    writer.close();

    return retval;
    // makeDalitzPlots(signalDalitz, "plots_from_mixfit/generated/");
}

GooPdf *makeBkg2_sigma() {
    // Variable* bkg2_sigma_js_meana = new Variable("bkg2_sigma_js_meana", 0.01, 0.01, -0.30, 0.30);
    // Variable* bkg2_sigma_js_sigma = new Variable("bkg2_sigma_js_sigma", 0.09, 0.01, 0, 0.4);
    // Variable* bkg2_sigma_js_gamma = new Variable("bkg2_sigma_js_gamma",-5.00, 0.10, -30, 0);
    // Variable* bkg2_sigma_js_delta = new Variable("bkg2_sigma_js_delta", 1.49, 0.01, 0.50, 5.00);
    // Variable* bkg2_sigma_frac_jsu = new Variable("bkg2_sigma_frac_jsu", 0.85, 0.01, 0.01, 1.00);
    // Variable* bkg2_sigma_frac_ga1 = new Variable("bkg2_sigma_frac_ga1", 0.04, 0.01, 0.01, 0.20);

    Variable bkg2_sigma_num_jsu("bkg2_sigma_num_jsu", 9100, 200, 1000, 15000);
    Variable bkg2_sigma_num_ga1("bkg2_sigma_num_ga1", 2400, 200, 500, 7000);
    Variable bkg2_sigma_num_ga2("bkg2_sigma_num_ga2", 2900, 200, 500, 7000);

    Variable bkg2_sigma_g1_meana("bkg2_sigma_g1_meana", 0.35, 0.01, 0.10, 0.50);
    Variable bkg2_sigma_g1_sigma("bkg2_sigma_g1_sigma", 0.30, 0.01, 0.05, 0.55);
    Variable bkg2_sigma_g2_meana("bkg2_sigma_g2_meana", 0.80, 0.01, 0.01, 1.50);
    Variable bkg2_sigma_g2_sigma("bkg2_sigma_g2_sigma", 0.90, 0.01, 0.01, 2.75);
    // JohnsonSUPdf* bkg2_sigma_js = new JohnsonSUPdf("bkg2_sigma_js", sigma, bkg2_sigma_js_meana, bkg2_sigma_js_sigma,
    // bkg2_sigma_js_gamma, bkg2_sigma_js_delta);

    Variable bkg2_sigma_js_meana("bkg2_sigma_js_meana", 0.35, 0.01, 0.00, 0.60);
    Variable bkg2_sigma_js_sigma("bkg2_sigma_js_sigma", 0.09, 0.01, 0, 0.4);
    Variable bkg2_sigma_js_gamma("bkg2_sigma_js_gamma", 2.00, 0.10, 0, 10);
    Variable bkg2_sigma_js_delta("bkg2_sigma_js_delta", 2);
    CrystalBallPdf *bkg2_sigma_js = new CrystalBallPdf(
        "bkg2_sigma_js", *sigma, bkg2_sigma_js_meana, bkg2_sigma_js_sigma, bkg2_sigma_js_gamma, bkg2_sigma_js_delta);

    GaussianPdf *bkg2_sigma_g1 = new GaussianPdf("bkg2_sigma_g1", *sigma, bkg2_sigma_g1_meana, bkg2_sigma_g1_sigma);
    GaussianPdf *bkg2_sigma_g2 = new GaussianPdf("bkg2_sigma_g2", *sigma, bkg2_sigma_g2_meana, bkg2_sigma_g2_sigma);

    weights.clear();
    weights.push_back(bkg2_sigma_num_jsu);
    weights.push_back(bkg2_sigma_num_ga1);
    weights.push_back(bkg2_sigma_num_ga2);
    comps.clear();
    comps.push_back(bkg2_sigma_js);
    comps.push_back(bkg2_sigma_g1);
    comps.push_back(bkg2_sigma_g2);

    GooPdf *ret = new AddPdf("bkg2_jsugg", weights, comps);
    return ret;
    // return bkg2_sigma_js;
}

GooPdf *makeBkg4_sigma() {
    // Variable* bkg4_sigma_js_meana = new Variable("bkg4_sigma_js_meana", 0.01, 0.01, -0.30, 0.30);
    // Variable* bkg4_sigma_js_sigma = new Variable("bkg4_sigma_js_sigma", 0.09, 0.01, 0, 0.4);
    // Variable* bkg4_sigma_js_gamma = new Variable("bkg4_sigma_js_gamma",-5.00, 0.10, -30, 0);
    // Variable* bkg4_sigma_js_delta = new Variable("bkg4_sigma_js_delta", 1.49, 0.01, 0.50, 5.00);
    // Variable* bkg4_sigma_frac_jsu = new Variable("bkg4_sigma_frac_jsu", 0.85, 0.01, 0.01, 1.00);
    // Variable* bkg4_sigma_frac_ga1 = new Variable("bkg4_sigma_frac_ga1", 0.04, 0.01, 0.01, 0.20);

    Variable bkg4_sigma_num_jsu("bkg4_sigma_num_jsu", 9100, 200, 1000, 15000);
    Variable bkg4_sigma_num_ga1("bkg4_sigma_num_ga1", 2400, 200, 500, 7000);
    Variable bkg4_sigma_num_ga2("bkg4_sigma_num_ga2", 2900, 200, 500, 7000);

    Variable bkg4_sigma_g1_meana("bkg4_sigma_g1_meana", 0.35, 0.01, 0.10, 0.50);
    Variable bkg4_sigma_g1_sigma("bkg4_sigma_g1_sigma", 0.30, 0.01, 0.05, 0.55);
    Variable bkg4_sigma_g2_meana("bkg4_sigma_g2_meana", 0.80, 0.01, 0.01, 1.50);
    Variable bkg4_sigma_g2_sigma("bkg4_sigma_g2_sigma", 0.90, 0.01, 0.01, 2.75);
    // JohnsonSUPdf* bkg4_sigma_js = new JohnsonSUPdf("bkg4_sigma_js", sigma, bkg4_sigma_js_meana, bkg4_sigma_js_sigma,
    // bkg4_sigma_js_gamma, bkg4_sigma_js_delta);

    Variable bkg4_sigma_js_meana("bkg4_sigma_js_meana", 0.35, 0.01, 0.00, 0.60);
    Variable bkg4_sigma_js_sigma("bkg4_sigma_js_sigma", 0.09, 0.01, 0, 0.4);
    Variable bkg4_sigma_js_gamma("bkg4_sigma_js_gamma", 2.00, 0.10, 0, 10);
    Variable bkg4_sigma_js_delta("bkg4_sigma_js_delta", 2);
    CrystalBallPdf *bkg4_sigma_js = new CrystalBallPdf(
        "bkg4_sigma_js", *sigma, bkg4_sigma_js_meana, bkg4_sigma_js_sigma, bkg4_sigma_js_gamma, bkg4_sigma_js_delta);

    GaussianPdf *bkg4_sigma_g1 = new GaussianPdf("bkg4_sigma_g1", *sigma, bkg4_sigma_g1_meana, bkg4_sigma_g1_sigma);
    GaussianPdf *bkg4_sigma_g2 = new GaussianPdf("bkg4_sigma_g2", *sigma, bkg4_sigma_g2_meana, bkg4_sigma_g2_sigma);

    weights.clear();
    weights.push_back(bkg4_sigma_num_jsu);
    weights.push_back(bkg4_sigma_num_ga1);
    weights.push_back(bkg4_sigma_num_ga2);
    comps.clear();
    comps.push_back(bkg4_sigma_js);
    comps.push_back(bkg4_sigma_g1);
    comps.push_back(bkg4_sigma_g2);

    GooPdf *ret = new AddPdf("bkg4_jsugg", weights, comps);
    return ret;
}

GooPdf *makeBkg3_sigma() {
    // Variable* bkg3_sigma_js_meana = new Variable("bkg3_sigma_js_meana", 0.05, 0.01, -0.30, 0.30);
    // Variable* bkg3_sigma_js_sigma = new Variable("bkg3_sigma_js_sigma", 0.013, 0.01, 0, 0.2);
    // Variable* bkg3_sigma_js_gamma = new Variable("bkg3_sigma_js_gamma",-6.00, 1.00, -30, 0);
    // Variable* bkg3_sigma_js_delta = new Variable("bkg3_sigma_js_delta", 1.99, 0.10, 0.50, 5.00);
    Variable bkg3_sigma_frac_jsu("bkg3_sigma_frac_jsu", 0.50, 0.01, 0.01, 1.00);
    Variable bkg3_sigma_frac_ga1("bkg3_sigma_frac_ga1", 0.04, 0.01, 0.01, 0.20);

    // Variable* bkg3_sigma_num_jsu = new Variable("bkg3_sigma_num_jsu", 11000, 200, 1000, 35000);
    // Variable* bkg3_sigma_num_ga1 = new Variable("bkg3_sigma_num_ga1", 9400, 200, 500, 10000);
    // Variable* bkg3_sigma_num_ga2 = new Variable("bkg3_sigma_num_ga2", 3900, 200, 500, 17000);
    // Variable* bkg3_sigma_num_ga3 = new Variable("bkg3_sigma_num_ga3", 3900, 200, 500, 7000);

    Variable bkg3_sigma_g1_meana("bkg3_sigma_g1_meana", 0.35, 0.01, 0.10, 0.50);
    Variable bkg3_sigma_g1_sigma("bkg3_sigma_g1_sigma", 0.10, 0.01, 0.01, 0.55);
    Variable bkg3_sigma_g2_meana("bkg3_sigma_g2_meana", 0.20, 0.01, 0.01, 1.50);
    Variable bkg3_sigma_g2_sigma("bkg3_sigma_g2_sigma", 0.10, 0.01, 0.001, 0.15);
    Variable bkg3_sigma_g2_gamma("bkg3_sigma_g2_gamma", -2.00, 1.00, -10, 10);
    Variable bkg3_sigma_g2_delta("bkg3_sigma_g2_delta", 2, 0.10, 0.50, 5.00);
    // Variable* bkg3_sigma_g3_meana = new Variable("bkg3_sigma_g3_meana", 0.20, 0.01, 0.01, 1.50);
    // Variable* bkg3_sigma_g3_sigma = new Variable("bkg3_sigma_g3_sigma", 0.20, 0.01, 0.01, 0.75);
    // Variable* bkg3_sigma_g2_width = new Variable("bkg3_sigma_g2_width", 0.10, 0.01, 0.01, 0.75);
    // JohnsonSUPdf* bkg3_sigma_js = new JohnsonSUPdf("bkg3_sigma_js", sigma, bkg3_sigma_js_meana, bkg3_sigma_js_sigma,
    // bkg3_sigma_js_gamma, bkg3_sigma_js_delta);

    Variable bkg3_sigma_js_meana("bkg3_sigma_js_meana", 0.35, 0.01, 0.00, 0.60);
    Variable bkg3_sigma_js_sigma("bkg3_sigma_js_sigma", 0.09, 0.01, 0, 0.40);
    Variable bkg3_sigma_js_gamma("bkg3_sigma_js_gamma", 2.00, 0.10, 0, 10);
    Variable bkg3_sigma_js_delta("bkg3_sigma_js_delta", 2);
    CrystalBallPdf *bkg3_sigma_js = new CrystalBallPdf(
        "bkg3_sigma_js", *sigma, bkg3_sigma_js_meana, bkg3_sigma_js_sigma, bkg3_sigma_js_gamma, bkg3_sigma_js_delta);
    // JohnsonSUPdf* bkg3_sigma_js = new JohnsonSUPdf("bkg3_sigma_js", sigma, bkg3_sigma_js_meana, bkg3_sigma_js_sigma,
    // bkg3_sigma_js_gamma, bkg3_sigma_js_delta);

    GaussianPdf *bkg3_sigma_g1 = new GaussianPdf("bkg3_sigma_g1", *sigma, bkg3_sigma_g1_meana, bkg3_sigma_g1_sigma);
    // GaussianPdf*  bkg3_sigma_g2 = new GaussianPdf("bkg3_sigma_g2", sigma, bkg3_sigma_g2_meana, bkg3_sigma_g2_sigma);
    // CrystalBallPdf* bkg3_sigma_g2 = new CrystalBallPdf("bkg3_sigma_g2", sigma, bkg3_sigma_g2_meana,
    // bkg3_sigma_g2_sigma, bkg3_sigma_g2_gamma, bkg3_sigma_g2_delta);
    JohnsonSUPdf *bkg3_sigma_g2 = new JohnsonSUPdf(
        "bkg3_sigma_g2", *sigma, bkg3_sigma_g2_meana, bkg3_sigma_g2_sigma, bkg3_sigma_g2_gamma, bkg3_sigma_g2_delta);
    // GaussianPdf*  bkg3_sigma_g3 = new GaussianPdf("bkg3_sigma_g3", sigma, bkg3_sigma_g3_meana, bkg3_sigma_g3_sigma);
    // VoigtianPdf*  bkg3_sigma_g2 = new VoigtianPdf("bkg3_sigma_g2", sigma, bkg3_sigma_g2_meana, bkg3_sigma_g2_sigma,
    // bkg3_sigma_g2_width);

    weights.clear();
    // weights.push_back(bkg3_sigma_num_jsu);
    // weights.push_back(bkg3_sigma_num_ga1);
    // weights.push_back(bkg3_sigma_num_ga2);
    // weights.push_back(bkg3_sigma_num_ga3);
    weights.push_back(bkg3_sigma_frac_jsu);
    weights.push_back(bkg3_sigma_frac_ga1);
    comps.clear();
    comps.push_back(bkg3_sigma_js);
    comps.push_back(bkg3_sigma_g1);
    comps.push_back(bkg3_sigma_g2);
    // comps.push_back(bkg3_sigma_g3);

    GooPdf *ret = new AddPdf("bkg3_jsugg", weights, comps);
    return ret;
}

GooPdf *makeGaussianTimePdf(int bkg) {
    // Gaussians for decay time.
    Variable *frac_ga2;
    Variable *frac_ga3;
    Variable *g1_meana;
    Variable *g1_sigma;
    Variable *g2_meana;
    Variable *g2_sigma;
    Variable *g3_meana;
    Variable *g3_sigma;

    std::string bkgname = "";

    switch(bkg) {
    case 4:
        frac_ga2 = new Variable("bkg4_frac_ga2", 0.60780, 0.01, 0.40, 0.80);
        frac_ga3 = new Variable("bkg4_frac_ga3", 0.04776, 0.01, 0.00, 0.20);
        g1_meana = new Variable("bkg4_g1_meana", 0.10164, 0.01, 0.00, 0.80);
        g1_sigma = new Variable("bkg4_g1_sigma", 0.27504, 0.01, 0.10, 0.80);
        g2_meana = new Variable("bkg4_g2_meana", 0.19974, 0.01, 0.10, 0.80);
        g2_sigma = new Variable("bkg4_g2_sigma", 0.63765, 0.01, 0.40, 0.80);
        g3_meana = new Variable("bkg4_g3_meana", 0.45817, 0.01, 0.20, 0.80);
        g3_sigma = new Variable("bkg4_g3_sigma", 1.52905, 0.01, 1.40, 1.80);
        bkgname  = "bkg4";
        break;

    case 3:
        frac_ga2 = new Variable("bkg3_frac_ga2", 0.51448, 0.01, 0.25, 0.75);
        frac_ga3 = new Variable("bkg3_frac_ga3", 0.04169, 0.01, 0.00, 0.40);
        g1_meana = new Variable("bkg3_g1_meana", 0.25101, 0.01, 0.01, 0.50);
        g1_sigma = new Variable("bkg3_g1_sigma", 0.31953, 0.01, 0.02, 0.50);
        g2_meana = new Variable("bkg3_g2_meana", 0.49139, 0.01, 0.25, 0.75);
        g2_sigma = new Variable("bkg3_g2_sigma", 0.65443, 0.01, 0.10, 1.00);
        g3_meana = new Variable("bkg3_g3_meana", 0.83600, 0.01, 0.50, 1.00);
        g3_sigma = new Variable("bkg3_g3_sigma", 1.51839, 0.01, 0.10, 2.00);
        bkgname  = "bkg3";
        break;

    case 2:
    default:
        frac_ga2 = new Variable("frac_ga2", 0.48994, 0.01, 0.1, 0.6);
        frac_ga3 = new Variable("frac_ga3", 0.04721, 0.01, 0.0, 0.1);
        g1_meana = new Variable("g1_meana", 0.01216, 0.01, -0.1, 0.1);
        g1_sigma = new Variable("g1_sigma", 0.25813, 0.01, 0.15, 0.35);
        g2_meana = new Variable("g2_meana", 0.05335, 0.01, -0.1, 0.1);
        g2_sigma = new Variable("g2_sigma", 0.58651, 0.01, 0.5, 1.2);
        g3_meana = new Variable("g3_meana", 0.17451, 0.01, 0.1, 1.85);
        g3_sigma = new Variable("g3_sigma", 1.15125, 0.01, 0.5, 1.3);
        bkgname  = "bkg2";
        break;
    }

    GaussianPdf *g1 = new GaussianPdf((bkgname + "_g1").c_str(), *dtime, *g1_meana, *g1_sigma);
    GaussianPdf *g2 = new GaussianPdf((bkgname + "_g2").c_str(), *dtime, *g2_meana, *g2_sigma);
    GaussianPdf *g3 = new GaussianPdf((bkgname + "_g3").c_str(), *dtime, *g3_meana, *g3_sigma);

    weights.clear();
    weights.push_back(*frac_ga2);

    if(3 == bkg)
        weights.push_back(*frac_ga3);

    comps.clear();
    comps.push_back(g2);

    if(3 == bkg)
        comps.push_back(g3);

    comps.push_back(g1);
    AddPdf *bkg_dtime = new AddPdf((bkgname + "_dtime").c_str(), weights, comps);
    return bkg_dtime;
}

GooPdf *makeExpGausTimePdf(int bkg) {
    std::string bkgname = "";

    switch(bkg) {
    case 4:
        bkgname = "bkg4";
        break;

    case 3:
        bkgname = "bkg3";
        break;

    case 2:
    default:
        bkgname = "bkg2";
        break;
    }

    Variable g1_mean((bkgname + "_dtime_gmean1"), 0, 0.01, -0.5, 0.5);
    Variable g1_sigm((bkgname + "_dtime_gsigm1"), 0.2, 0.01, 0.01, 0.8);
    Variable e1_alph((bkgname + "_dtime_alpha1"), 2.5, 0.01, 0.01, 7.5);

    Variable g2_mean((bkgname + "_dtime_gmean2"), -0.3, 0.01, -0.85, 0.85);
    Variable g2_sigm((bkgname + "_dtime_gsigm2"), 0.2, 0.01, 0.01, 0.8);
    Variable e2_alph((bkgname + "_dtime_alpha2"), 0.5, 0.01, 0.01, 10.0);

    ExpGausPdf *exp1 = new ExpGausPdf((bkgname + "_exp1"), *dtime, g1_mean, g1_sigm, e1_alph);
    ExpGausPdf *exp2 = new ExpGausPdf((bkgname + "_exp2"), *dtime, g2_mean, g2_sigm, e2_alph);

    Variable frac1((bkgname + "_dtime_frac1"), 0.1, 0.01, 0, 0.8);

    GooPdf *ret = new AddPdf((bkgname + "_dtime"), frac1, exp1, exp2);

    return ret;
}

GooPdf *makeBkg2DalitzPdf(bool fixem = true) {
    if(!kzero_veto)
        makeKzeroVeto();

    GooPdf *bkg2_dalitz = nullptr;

    if(Parameter == bkg2Model) {
        comps.clear();

        vector<Variable> offsets;
        vector<Observable> observables;
        vector<Variable> coefficients;
        offsets.push_back(constantOne);
        offsets.push_back(constantOne);
        observables.push_back(*m12);
        observables.push_back(*m13);
        double weightOffset = 3;
        // Recurring factor 3 offsets division by total weight in AddPdf.
        coefficients.push_back(Variable("bkg2_x0y0", 1.0 * weightOffset));
        coefficients.push_back(
            Variable("bkg2_x1y0", 0.13184 * weightOffset, 0.01, 0.01 * weightOffset, 0.18 * weightOffset));
        coefficients.push_back(
            Variable("bkg2_x2y0", 0.02062 * weightOffset, 0.01, 0.00 * weightOffset, 0.17 * weightOffset));
        coefficients.push_back(
            Variable("bkg2_x3y0", 0.04688 * weightOffset, 0.01, 0.00 * weightOffset, 0.08 * weightOffset));
        coefficients.push_back(
            Variable("bkg2_x0y1", -0.02568 * weightOffset, 0.01, -0.15 * weightOffset, 0.04 * weightOffset));
        coefficients.push_back(
            Variable("bkg2_x1y1", 0.06805 * weightOffset, 0.01, 0.02 * weightOffset, 0.10 * weightOffset));
        coefficients.push_back(
            Variable("bkg2_x2y1", 0.38557 * weightOffset, 0.01, 0.30 * weightOffset, 0.50 * weightOffset));
        coefficients.push_back(
            Variable("bkg2_x0y2", 0.11252 * weightOffset, 0.01, 0.05 * weightOffset, 0.20 * weightOffset));
        coefficients.push_back(
            Variable("bkg2_x1y2", 0.24896 * weightOffset, 0.01, 0.20 * weightOffset, 0.50 * weightOffset));
        coefficients.push_back(
            Variable("bkg2_x0y3", 0.05605 * weightOffset, 0.01, -0.05 * weightOffset, 0.15 * weightOffset));

        PolynomialPdf *poly = new PolynomialPdf("bkg2Pdf", observables, coefficients, offsets, 3);

        Variable bkg2_decZmin("bkg2_decZmin", 3.30411);
        Variable bkg2_conZmin("bkg2_conZmin", 0.29909);
        TrigThresholdPdf *bkg2_loZ
            = new TrigThresholdPdf("bkg2_loZ", *m12, *m13, minDalitzZ, bkg2_decZmin, bkg2_conZmin, massSum, false);
        // bkg2_loZ->setDebugMask(1);

        comps.clear();
        comps.push_back(poly);
        comps.push_back(bkg2_loZ);
        comps.push_back(kzero_veto);
        // Separate PDF to avoid triggering numerical normalisation over all four observables.
        ProdPdf *poly_x_veto = new ProdPdf("poly_x_veto", comps);
        // poly_x_veto->setDebugMask(1);

        // One omega->pipipi0 reflection.
        // Factor 3 in amplitudes is to offset division by total weight in AddPdf.
        DecayInfo3 special_rho_decay;
        special_rho_decay.motherMass   = _mD0;
        special_rho_decay.daug1Mass    = piZeroMass;
        special_rho_decay.daug2Mass    = piPlusMass;
        special_rho_decay.daug3Mass    = piPlusMass;
        special_rho_decay.meson_radius = 1.5;

        ResonancePdf *bkg2_rho_ref = new Resonances::Gauss(
            "bkg2_rho_ref",
            Variable("bkg2_rho_ref_amp", 0.00896 * weightOffset, 0.001, 0, 0.015 * weightOffset),
            constantZero,
            Variable("bkg2_rho_ref_mass", 0.53172),
            Variable("bkg2_rho_ref_width", 0.06426),
            PAIR_13);
        special_rho_decay.resonances.push_back(bkg2_rho_ref);

        Variable bkg2_rho_poly_offset("bkg2_rho_poly_offset", 1.64254);
        Variable bkg2_rho_poly_linear("bkg2_rho_poly_linear", 0);
        Variable bkg2_rho_poly_second("bkg2_rho_poly_second", -0.48166);
        weights.clear();
        weights.push_back(constantOne);
        weights.push_back(bkg2_rho_poly_linear);
        weights.push_back(bkg2_rho_poly_second);
        PolynomialPdf *bkg2_rho_poly = new PolynomialPdf("bkg2_rho_poly", *m12, weights, bkg2_rho_poly_offset);
        comps.clear();
        comps.push_back(kzero_veto);
        comps.push_back(bkg2_rho_poly);
        comps.push_back(bkg2_loZ);
        ProdPdf *bkg2_rho_mods = new ProdPdf("bkg2_rho_mods", comps);
        incsum1 = new IncoherentSumPdf("incsum1", *m12, *m13, *eventNumber, special_rho_decay, bkg2_rho_mods);

        // Three spin-0 rho resonances to be added incoherently.
        DecayInfo3 incoherent_rho0s;
        incoherent_rho0s.motherMass   = _mD0;
        incoherent_rho0s.daug1Mass    = piZeroMass;
        incoherent_rho0s.daug2Mass    = piPlusMass;
        incoherent_rho0s.daug3Mass    = piPlusMass;
        incoherent_rho0s.meson_radius = 0; // Mikhail uses zero radius for incoherent resonances.

        ResonancePdf *bkg2_incRho0 = new Resonances::RBW(
            "bkg2_incRho0",
            Variable("bkg2_incRho0_amp", 0.00304 * weightOffset, 0.001, 0.0, 0.006 * weightOffset),
            constantZero,
            fixedRhoMass,
            fixedRhoWidth,
            0, // Incoherent rho has effective spin 0.
            PAIR_23);
        incoherent_rho0s.resonances.push_back(bkg2_incRho0);

        ResonancePdf *bkg2_incRhoP = new Resonances::RBW(
            "bkg2_incRhoP",
            Variable("bkg2_incRhoP_amp", 0.00586 * weightOffset, 0.001, 0.0, 0.012 * weightOffset),
            constantZero,
            fixedRhoMass,
            fixedRhoWidth,
            0,
            PAIR_12);
        incoherent_rho0s.resonances.push_back(bkg2_incRhoP);

        ResonancePdf *bkg2_incRhoM = new Resonances::RBW(
            "bkg2_incRhoM",
            Variable("bkg2_incRhoM_amp", 0.00635 * weightOffset, 0.001, 0.0, 0.015 * weightOffset),
            constantZero,
            fixedRhoMass,
            fixedRhoWidth,
            0,
            PAIR_13);
        incoherent_rho0s.resonances.push_back(bkg2_incRhoM);

        comps.clear();
        comps.push_back(kzero_veto);
        comps.push_back(bkg2_loZ);
        ProdPdf *bkg2_rho_mods2 = new ProdPdf("bkg2_rho_mods2", comps);

        incsum2 = new IncoherentSumPdf("incsum2", *m12, *m13, *eventNumber, incoherent_rho0s, bkg2_rho_mods2);

        weights.clear();
        weights.push_back(constantOne);
        weights.push_back(constantOne);
        weights.push_back(constantOne);
        comps.clear();
        comps.push_back(poly_x_veto);
        comps.push_back(incsum1);
        comps.push_back(incsum2);

        bkg2_dalitz = new AddPdf("bkg2_dalitz", weights, comps);
        bkg2_dalitz->addSpecialMask(PdfBase::ForceCommonNorm);
        // bkg2_dalitz->setDebugMask(1);

    } else if(Histogram == bkg2Model) {
        bkg2_dalitz = makeBackgroundHistogram(2);
    } else if(Sideband == bkg2Model) {
        comps.clear();
        // Originally had massD0 sidebands separated into deltaM high and low,
        // but these distributions were extremely similar - collapsed them
        // into just massD0 sidebands.
        std::string fname = app_ptr->get_filename("dataFiles/sideband1.txt", "examples/pipipi0DPFit");
        comps.push_back(makeBackgroundHistogram(101, fname));
        fname = app_ptr->get_filename("dataFiles/sideband2.txt", "examples/pipipi0DPFit");
        comps.push_back(makeBackgroundHistogram(102, fname));
        // comps.push_back(makeBackgroundHistogram(103, "./dataFiles/sideband3.txt"));
        // comps.push_back(makeBackgroundHistogram(104, "./dataFiles/sideband4.txt"));
        weights.clear();
        weights.push_back(Variable("sband1Weight", 300000, 1000, 100, 750000));
        weights.push_back(Variable("sband2Weight", 100000, 1000, 100, 500000));
        // weights.push_back(new Variable("sband3Weight", 150000, 1000, 100, 500000));
        // weights.push_back(new Variable("sband4Weight",  50000, 1000, 100, 500000));
        bkg2_dalitz = new AddPdf("bkg2_dalitz", weights, comps);
    } else {
        // This cannot possibly happen, and if it does something is wrong. Panic!
        assert(Sideband == bkg2Model);
    }

    GooPdf *bkg2_dtime = nullptr;

    if(gaussBkgTime)
        bkg2_dtime = makeGaussianTimePdf(2);
    else
        bkg2_dtime = makeExpGausTimePdf(2);

    // bkg2_jsugg = sig0_jsugg; // Mikhail uses same sigma distribution for everything.
    // Separate sigma_t distribution
    // bkg2_jsugg = makeBkg2_sigma();
    bkg2_jsugg = makeBkg_sigma_strips(2);
    bkg2_jsugg->addSpecialMask(PdfBase::ForceSeparateNorm);

    // Finally create overall product.
    comps.clear();
    // poly->setDebugMask(1);
    // bkg2_dalitz->setDebugMask(1);
    // incsum2->setDebugMask(1);
    comps.push_back(bkg2_dalitz);
    comps.push_back(bkg2_dtime);
    comps.push_back(bkg2_jsugg);

    GooPdf *ret = new ProdPdf("bkg2_total", comps);

    if(fixem)
        ret->setParameterConstantness(true);

    return ret;
}

GooPdf *makeBkg3Eff() {
    // Smoothed histogram from flat-file data.
    // Only 4500 events, so use large bins.

    obsweights.clear();
    obsweights.push_back(*m12);
    obsweights.push_back(*m13);

    int m12bins = m12->getNumBins();
    int m13bins = m13->getNumBins();

    m12->setNumBins(30);
    m13->setNumBins(30);
    BinnedDataSet *bkg3_eff_data = new BinnedDataSet(obsweights);
    std::ifstream reader;
    std::string fname = app_ptr->get_filename("dataFiles/efficiency_bkg3_flat.txt", "examples/pipipi0DPFit");
    readWrapper(reader, fname);
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

        if(reader.eof())
            break;

        reader >> dummy; // m23, m(pi+ pi-), called m12 in processToyRoot convention.
        reader >> *m12;  // Already swapped according to D* charge
        reader >> *m13;

        // Everything else is irrelevant for this purpose!
        for(int i = 0; i < 16; ++i)
            reader >> dummy;

        bkg3_eff_data->addEvent();
    }

    Variable bkg3_eff_smoothing("bkg3_eff_smoothing", 1.0, 0.1, 0, 1.25);
    // Variable* bkg3_eff_smoothing = new Variable("bkg3_eff_smoothing", 1.0);
    SmoothHistogramPdf *ret = new SmoothHistogramPdf("bkg3_efficiency", bkg3_eff_data, bkg3_eff_smoothing);

    m12->setNumBins(m12bins);
    m13->setNumBins(m13bins);

    return ret;
}

SmoothHistogramPdf *makeBackgroundHistogram(int bkgnum, std::string overridename) {
    std::ifstream reader;
    sprintf(strbuffer, "./dataFiles/bkgDalitz_%i.txt", bkgnum);

    if(overridename != "")
        sprintf(strbuffer, "%s", overridename.c_str());

    std::string fname = app_ptr->get_filename(strbuffer, "examples/pipipi0DPFit");
    readWrapper(reader, fname);
    std::string buffer;

    while(!reader.eof()) {
        reader >> buffer;

        if(buffer == "====")
            break;

        std::cout << buffer;
    }

    obsweights.clear();
    obsweights.push_back(*m12);
    obsweights.push_back(*m13);
    BinnedDataSet *bkg_binned_data = new BinnedDataSet(obsweights);

    double dummy = 0;

    while(!reader.eof()) {
        reader >> dummy;
        reader >> dummy; // m23, m(pi+ pi-), called m12 in processToyRoot convention.
        reader >> *m12;  // Already swapped according to D* charge. m12 = m(pi+pi0)
        reader >> *m13;

        // Don't need the rest.
        for(int i = 0; i < 16; ++i)
            reader >> dummy;

        bkg_binned_data->addEvent();

        // std::cout << m12->getValue() << " " << m13->getValue() << std::endl;
    }

    std::cout << "Read " << bkg_binned_data->getNumEvents() << " events for background " << bkgnum << std::endl;
    sprintf(strbuffer, "bkg%i_dalitz_smoothing", bkgnum);
    Variable smoothing(strbuffer, 1);

    if((-1 != bkgHistRandSeed) && ((3 == bkgnum) || (4 == bkgnum))) {
        std::cout << "Shuffling background " << bkgnum << " histogram with random seed " << bkgHistRandSeed
                  << std::endl;
        TRandom donram(bkgHistRandSeed);

        for(unsigned int bin = 0; bin < bkg_binned_data->getNumBins(); ++bin) {
            double events = bkg_binned_data->getBinContent(bin);

            if(1 > events)
                continue;

            double newEvents = -1;

            while(0 > newEvents)
                newEvents = donram.Gaus(events, sqrt(events));

            bkg_binned_data->setBinContent(bin, newEvents);
        }
    }

    sprintf(strbuffer, "bkg%i_dalitz", bkgnum);
    SmoothHistogramPdf *bkg_dalitz = new SmoothHistogramPdf(strbuffer, bkg_binned_data, smoothing);
    // bkg_dalitz->setDebugMask(1);
    return bkg_dalitz;
}

GooPdf *makeBackground3DalitzParam() {
    // I can't make this thing describe the background 3 data.

    // GooPdf* bkg3_eff = makeBkg3Eff();
    weights.clear();

    vector<Variable> offsets;
    vector<Observable> observables;
    vector<Variable> coefficients;
    offsets.push_back(constantOne);
    offsets.push_back(constantOne);

    observables.push_back(*m12);
    observables.push_back(*m13);
    // Recurring factor 3 offsets division by total weight in AddPdf.
    double weightOffset = 1;

    coefficients.push_back(Variable("bkg3_x0y0", 1.00 * weightOffset));
    // coefficients.push_back(new Variable("bkg3_x0y0",  1.10 * weightOffset, 0.01, 0.01 * weightOffset, 1.50 *
    // weightOffset));
    coefficients.push_back(
        Variable("bkg3_x1y0", -0.36937 * weightOffset, 0.01, -1.50 * weightOffset, 0.00 * weightOffset));
    coefficients.push_back(
        Variable("bkg3_x2y0", 1.36184 * weightOffset, 0.01, -0.10 * weightOffset, 1.60 * weightOffset));
    // coefficients.push_back(new Variable("bkg3_x3y0", -0.43177 * weightOffset, 0.01,-1.60*weightOffset,
    // 0.60*weightOffset));
    coefficients.push_back(
        Variable("bkg3_x0y1", -0.27691 * weightOffset, 0.01, -1.50 * weightOffset, 0.00 * weightOffset));
    coefficients.push_back(
        Variable("bkg3_x1y1", 2.16029 * weightOffset, 0.01, 0.30 * weightOffset, 4.50 * weightOffset));
    // coefficients.push_back(new Variable("bkg3_x2y1", -2.04133 * weightOffset, 0.01,-2.50*weightOffset,
    // 1.50*weightOffset));
    coefficients.push_back(
        Variable("bkg3_x0y2", 1.33100 * weightOffset, 0.01, 1.00 * weightOffset, 2.00 * weightOffset));
    // coefficients.push_back(new Variable("bkg3_x1y2", -1.88226 * weightOffset, 0.01,-2.20*weightOffset,
    // 1.00*weightOffset));
    // coefficients.push_back(new Variable("bkg3_x0y3", -0.58920 * weightOffset, 0.01,-1.00*weightOffset,
    // 2.00*weightOffset));
    // PolynomialPdf* poly = new PolynomialPdf("bkg3Pdf", observables, coefficients, offsets, 3);
    PolynomialPdf *poly = new PolynomialPdf("bkg3Pdf", observables, coefficients, offsets, 2);

    // coefficients.push_back(new Variable("bkg3_x0y0",  0.10 * weightOffset, 0.01, 0.01 * weightOffset, 0.20 *
    // weightOffset));
    // PolynomialPdf* poly = new PolynomialPdf("bkg3Pdf", observables, coefficients, offsets, 0);

    // Background 3 does not have a trig threshold in Mikhail's fit
    // - that is, it has one, but the dec variable is set above the
    // threshold, so it doesn't do anything.
    // That's Mikhail's fit; I'm putting one in to try to deal with
    // the asymmetry in the rho+.
    // Didn't work.
    // Variable* bkg3_decZmin = new Variable("bkg3_decZmin", 3.30411, 0.1, 1, 5);
    // Variable* bkg3_conZmin = new Variable("bkg3_conZmin", 0.29909, 0.01, 0.1, 0.9);
    // TrigThresholdPdf* bkg3_loZ = new TrigThresholdPdf("bkg3_loZ", m12, m13, minDalitzZ, bkg3_decZmin, bkg3_conZmin,
    // massSum, false);

    comps.clear();
    comps.push_back(poly);
    comps.push_back(kzero_veto);
    // comps.push_back(bkg3_eff);
    // comps.push_back(bkg3_loZ);

    ProdPdf *poly_x_veto = new ProdPdf("poly_x_veto", comps);

    // One misIDpi0.
    // Factor 3 in amplitudes is to offset division by total weight in AddPdf.
    DecayInfo3 special_pi0_decay;
    special_pi0_decay.motherMass   = _mD0;
    special_pi0_decay.daug1Mass    = piZeroMass;
    special_pi0_decay.daug2Mass    = piPlusMass;
    special_pi0_decay.daug3Mass    = piPlusMass;
    special_pi0_decay.meson_radius = 1.5;

    ResonancePdf *bkg3_pi0_ref = new Resonances::Gauss(
        "bkg3_pi0_ref",
        Variable("bkg3_pi0_ref_amp", 0.01189 * weightOffset, 0.01, 0.00 * weightOffset, 0.25 * weightOffset),
        constantZero,
        Variable("bkg3_pi0_ref_mass", 1.65766, 0.01, 1.4, 1.8),
        Variable("bkg3_pi0_ref_width", 0.05018, 0.01, 0.02, 0.20),
        PAIR_23);
    special_pi0_decay.resonances.push_back(bkg3_pi0_ref);

    // Mikhail defines 'transverse Z' as y - x - (parameter).
    // Variable* bkg3_pi0_transZ_offset = new Variable("bkg3_pi0_transZ_offset", -0.04381, 0.001, -0.5, 0.5);
    Variable bkg3_pi0_transZ_offset("bkg3_pi0_transZ_offset", -0.04381);
    offsets.clear();
    observables.clear();
    coefficients.clear();
    observables.push_back(*m12);
    observables.push_back(*m13);
    coefficients.push_back(bkg3_pi0_transZ_offset);
    coefficients.push_back(constantMinusOne);
    coefficients.push_back(constantOne);
    offsets.push_back(constantZero);
    offsets.push_back(constantZero);
    PolynomialPdf *bkg3_pi0_transZ = new PolynomialPdf("bkg3_pi0_transZ", observables, coefficients, offsets, 1);

    // Now we're going to take (1 - tz^2 * (parameter)) and multiply that into the misID pi0.
    // Variable* bkg3_pi0_transZ_quad = new Variable("bkg3_pi0_transZ_quad", 2.12277, 0.01, -1.5, 6.0);
    Variable bkg3_pi0_transZ_quad("bkg3_pi0_transZ_quad", 2.12277);
    coefficients.clear();
    coefficients.push_back(constantOne);
    coefficients.push_back(constantZero);
    coefficients.push_back(bkg3_pi0_transZ_quad);
    // Notice the fake dependence of the polynomial on m12; in fact CompositePdf
    // will send it a fake event, we just have to supply a reasonable index.
    PolynomialPdf *bkg3_pi0_shell       = new PolynomialPdf("bkg3_pi0_shell", *m12, coefficients);
    CompositePdf *bkg3_pi0_transZ_total = new CompositePdf("bkg3_pi0_transZ_total", bkg3_pi0_transZ, bkg3_pi0_shell);

    comps.clear();
    comps.push_back(kzero_veto);
    comps.push_back(bkg3_pi0_transZ_total);
    // comps.push_back(bkg3_eff);
    // comps.push_back(bkg3_loZ);
    // ProdPdf* bkg3_pi0_mods = new ProdPdf("bkg3_pi0_mods", comps);
    // incsum3 = new IncoherentSumPdf("incsum3", m12, m13, eventNumber, special_pi0_decay, bkg3_pi0_mods);

    // Three spin-1 rho resonances to be added incoherently.
    DecayInfo3 incoherent_rhos;
    incoherent_rhos.motherMass   = _mD0;
    incoherent_rhos.daug1Mass    = piZeroMass;
    incoherent_rhos.daug2Mass    = piPlusMass;
    incoherent_rhos.daug3Mass    = piPlusMass;
    incoherent_rhos.meson_radius = 0; // Mikhail uses zero radius for incoherent resonances.

    ResonancePdf *bkg3_incRho0 = new Resonances::RBW(
        "bkg3_incRho0",
        Variable("bkg3_incRho0_amp", 0.00807 * weightOffset, 0.01, 0.00 * weightOffset, 0.25 * weightOffset),
        constantZero,
        Variable("bkg3_incRho0_mass", 0.800, 0.01, 0.6, 1.0),
        Variable("bkg3_incRho0_width", 0.15, 0.01, 0.10, 0.40),
        1, // These rhos are spin 1, being bad signal.
        PAIR_23);
    incoherent_rhos.resonances.push_back(bkg3_incRho0);

    ResonancePdf *bkg3_incRhoP = new Resonances::RBW(
        "bkg3_incRhoP",
        Variable("bkg3_incRhoP_amp", 0.01683 * weightOffset, 0.01, 0.00 * weightOffset, 0.25 * weightOffset),
        constantZero,
        Variable("bkg3_incRhoP_mass", 0.800, 0.01, 0.6, 1.0),
        Variable("bkg3_incRhoP_width", 0.15, 0.01, 0.10, 0.40),
        1,
        PAIR_12);
    incoherent_rhos.resonances.push_back(bkg3_incRhoP);

    ResonancePdf *bkg3_incRhoM = new Resonances::RBW(
        "bkg3_incRhoM",
        Variable("bkg3_incRhoM_amp", 0.01645 * weightOffset, 0.01, 0.00 * weightOffset, 0.25 * weightOffset),
        constantZero,
        Variable("bkg3_incRhoM_mass", 0.900, 0.01, 0.6, 1.0),
        Variable("bkg3_incRhoM_width", 0.35, 0.01, 0.10, 0.60),
        1,
        PAIR_13);
    incoherent_rhos.resonances.push_back(bkg3_incRhoM);

    comps.clear();
    comps.push_back(kzero_veto);
    // comps.push_back(bkg3_loZ);
    // comps.push_back(bkg3_eff);
    // ProdPdf* bkg3_rho_mods = new ProdPdf("bkg3_rho_mods", comps);

    // incsum4 = new IncoherentSumPdf("incsum4", m12, m13, eventNumber, incoherent_rhos, bkg3_rho_mods);
    // incsum4 = new IncoherentSumPdf("incsum4", m12, m13, eventNumber, incoherent_rhos, kzero_veto);

    weights.clear();
    weights.push_back(constantOne);
    // weights.push_back(constantOne);
    // weights.push_back(constantOne);
    comps.clear();
    comps.push_back(poly_x_veto);
    // comps.push_back(incsum3);
    // comps.push_back(incsum4);

    AddPdf *bkg3_dalitz = new AddPdf("bkg3_dalitz", weights, comps);
    bkg3_dalitz->addSpecialMask(PdfBase::ForceCommonNorm);
    return bkg3_dalitz;
}

GooPdf *makeBackground4DalitzParam() {
    vector<Variable> offsets;
    vector<Observable> observables;
    vector<Variable> coefficients;
    offsets.push_back(constantOne);
    offsets.push_back(constantOne);

    observables.push_back(*m12);
    observables.push_back(*m13);
    // Recurring factor 3 offsets division by total weight in AddPdf.

    double weightOffset = 3;
    // coefficients.push_back(new Variable("bkg4_x0y0",  1.0 * weightOffset, 0.01, 0.50*weightOffset,
    // 1.50*weightOffset));
    coefficients.push_back(Variable("bkg4_x0y0", 1.0 * weightOffset));
    coefficients.push_back(
        Variable("bkg4_x1y0", -0.18594 * weightOffset, 0.01, -0.50 * weightOffset, 0.50 * weightOffset));
    coefficients.push_back(
        Variable("bkg4_x2y0", 0.45459 * weightOffset, 0.01, 0.25 * weightOffset, 0.75 * weightOffset));
    coefficients.push_back(
        Variable("bkg4_x3y0", -0.20869 * weightOffset, 0.01, -0.50 * weightOffset, 0.50 * weightOffset));
    coefficients.push_back(
        Variable("bkg4_x0y1", -0.65061 * weightOffset, 0.01, -1.50 * weightOffset, 0.50 * weightOffset));
    coefficients.push_back(
        Variable("bkg4_x1y1", 0.11000 * weightOffset, 0.01, 0.00 * weightOffset, 0.50 * weightOffset));
    coefficients.push_back(
        Variable("bkg4_x2y1", 0.42009 * weightOffset, 0.01, 0.25 * weightOffset, 1.00 * weightOffset));
    coefficients.push_back(
        Variable("bkg4_x0y2", -0.06151 * weightOffset, 0.01, -0.50 * weightOffset, 0.50 * weightOffset));
    coefficients.push_back(
        Variable("bkg4_x1y2", 0.58508 * weightOffset, 0.01, 0.20 * weightOffset, 1.50 * weightOffset));
    coefficients.push_back(
        Variable("bkg4_x0y3", 0.54740 * weightOffset, 0.01, 0.20 * weightOffset, 1.50 * weightOffset));

    PolynomialPdf *poly = new PolynomialPdf("bkg4Pdf", observables, coefficients, offsets, 3);

    Variable bkg4_decZmin("bkg4_decZmin", 2.77576);
    Variable bkg4_conZmin("bkg4_conZmin", 0.23328);
    TrigThresholdPdf *bkg4_loZ
        = new TrigThresholdPdf("bkg4_loZ", *m12, *m13, minDalitzZ, bkg4_decZmin, bkg4_conZmin, massSum, false);

    comps.clear();
    comps.push_back(poly);
    comps.push_back(bkg4_loZ);
    comps.push_back(kzero_veto);
    // Separate PDF to avoid triggering numerical normalisation over all four observables.
    ProdPdf *poly_x_veto = new ProdPdf("poly_x_veto", comps);

    // One pipi bump.
    // Factor 3 in amplitudes is to offset division by total weight in AddPdf.
    DecayInfo3 special_pipi_decay;
    special_pipi_decay.motherMass   = _mD0;
    special_pipi_decay.daug1Mass    = piZeroMass;
    special_pipi_decay.daug2Mass    = piPlusMass;
    special_pipi_decay.daug3Mass    = piPlusMass;
    special_pipi_decay.meson_radius = 1.5;

    ResonancePdf *bkg4_pipi_ref = new Resonances::Gauss("bkg4_pipi_ref",
                                                        Variable("bkg4_pipi_ref_amp", 0.00147 * weightOffset),
                                                        constantZero,
                                                        Variable("bkg4_pipi_ref_mass", 1.32447),
                                                        Variable("bkg4_pipi_ref_width", 0.04675),
                                                        PAIR_23);
    special_pipi_decay.resonances.push_back(bkg4_pipi_ref);

    // Mikhail defines 'transverse Z' as y - x - (parameter).
    Variable bkg4_pipi_transZ_offset("bkg4_pipi_transZ_offset", -0.39877);

    offsets.clear();
    observables.clear();
    coefficients.clear();
    observables.push_back(*m12);
    observables.push_back(*m13);
    coefficients.push_back(bkg4_pipi_transZ_offset);
    coefficients.push_back(constantMinusOne);
    coefficients.push_back(constantOne);
    offsets.push_back(constantZero);
    offsets.push_back(constantZero);
    PolynomialPdf *bkg4_pipi_transZ = new PolynomialPdf("bkg4_pipi_transZ", observables, coefficients, offsets, 1);

    // Now we're going to take (1 - tz^2 * (parameter)) and multiply that into the pipi bump.
    Variable bkg4_pipi_transZ_quad("bkg4_pipi_transZ_quad", -0.25640);
    coefficients.clear();
    coefficients.push_back(constantOne);
    coefficients.push_back(constantZero);
    coefficients.push_back(bkg4_pipi_transZ_quad);
    // Notice the fake dependence of the polynomial on m12; in fact CompositePdf
    // will send it a fake event, we just have to supply a reasonable index.
    PolynomialPdf *bkg4_pipi_shell = new PolynomialPdf("bkg4_pipi_shell", *m12, coefficients);
    CompositePdf *bkg4_pipi_transZ_total
        = new CompositePdf("bkg4_pipi_transZ_total", bkg4_pipi_transZ, bkg4_pipi_shell);

    comps.clear();
    comps.push_back(kzero_veto);
    comps.push_back(bkg4_loZ);
    comps.push_back(bkg4_pipi_transZ_total);

    ProdPdf *bkg4_pipi_mods = new ProdPdf("bkg4_pipi_mods", comps);
    incsum5 = new IncoherentSumPdf("incsum5", *m12, *m13, *eventNumber, special_pipi_decay, bkg4_pipi_mods);

    // Three spin-0 rho resonances to be added incoherently.
    DecayInfo3 incoherent_rho0s;
    incoherent_rho0s.motherMass   = _mD0;
    incoherent_rho0s.daug1Mass    = piZeroMass;
    incoherent_rho0s.daug2Mass    = piPlusMass;
    incoherent_rho0s.daug3Mass    = piPlusMass;
    incoherent_rho0s.meson_radius = 0; // Mikhail uses zero radius for incoherent resonances.

    ResonancePdf *bkg4_incRho0 = new Resonances::RBW("bkg4_incRho0",
                                                     Variable("bkg4_incRho0_amp", 0.00429 * weightOffset),
                                                     constantZero,
                                                     fixedRhoMass,
                                                     fixedRhoWidth,
                                                     0, // These rhos are spin 0.
                                                     PAIR_23);
    incoherent_rho0s.resonances.push_back(bkg4_incRho0);

    ResonancePdf *bkg4_incRhoP = new Resonances::RBW("bkg4_incRhoP",
                                                     Variable("bkg4_incRhoP_amp", 0.00705 * weightOffset),
                                                     constantZero,
                                                     fixedRhoMass,
                                                     fixedRhoWidth,
                                                     0,
                                                     PAIR_12);
    incoherent_rho0s.resonances.push_back(bkg4_incRhoP);

    ResonancePdf *bkg4_incRhoM = new Resonances::RBW("bkg4_incRhoM",
                                                     Variable("bkg4_incRhoM_amp", -0.00043 * weightOffset),
                                                     constantZero,
                                                     fixedRhoMass,
                                                     fixedRhoWidth,
                                                     0,
                                                     PAIR_13);
    incoherent_rho0s.resonances.push_back(bkg4_incRhoM);

    comps.clear();
    comps.push_back(kzero_veto);
    comps.push_back(bkg4_loZ);
    ProdPdf *bkg4_incrho_mods = new ProdPdf("bkg4_incrho_mods", comps);
    incsum6 = new IncoherentSumPdf("incsum6", *m12, *m13, *eventNumber, incoherent_rho0s, bkg4_incrho_mods);
    // incsum6 = new IncoherentSumPdf("incsum6", m12, m13, eventNumber, incoherent_rho0s, kzero_veto);

    weights.clear();
    weights.push_back(constantOne);
    weights.push_back(constantOne);
    weights.push_back(constantOne);
    comps.clear();
    comps.push_back(poly_x_veto);
    comps.push_back(incsum5);
    comps.push_back(incsum6);

    AddPdf *bkg4_dalitz = new AddPdf("bkg4_dalitz", weights, comps);
    bkg4_dalitz->addSpecialMask(PdfBase::ForceCommonNorm);
    return bkg4_dalitz;
}

GooPdf *makeBkg3DalitzPdf(bool fixem = true) {
    if(!kzero_veto)
        makeKzeroVeto();

    comps.clear();
    weights.clear();

    // Using a histogram for Dalitz description. Notice: The purpose of this
    // fit is to get a description for use in fitting actual data. When fitting
    // data I can just use the background 3 Monte Carlo histogram. When fitting
    // the MC, it doesn't matter what I do, because I'm going to be using the
    // histogram. So I load up all the MC data and use it in the histogram, either
    // way.

    GooPdf *bkg3_dalitz = nullptr;

    if(!notUseBackground3Hist)
        bkg3_dalitz = makeBackgroundHistogram(3);
    else
        bkg3_dalitz = makeBackground3DalitzParam();

    // bkg3_dalitz->setDebugMask(1);

    GooPdf *bkg3_dtime = nullptr;

    if(gaussBkgTime)
        bkg3_dtime = makeGaussianTimePdf(3);
    else
        bkg3_dtime = makeExpGausTimePdf(3);

    // bkg3_jsugg = makeBkg3_sigma();
    // bkg3_jsugg = sig0_jsugg; // Mikhail uses same sigma distribution for everything.
    bkg3_jsugg = makeBkg_sigma_strips(3);
    bkg3_jsugg->addSpecialMask(PdfBase::ForceSeparateNorm);
    // Otherwise ProdPdf tries to use the default overall integration,
    // because bkg3_jsugg depends on m12, m13 due to the striping, and that has
    // disastrous results for bkg3_dalitz. Note that this doesn't, actually,
    // contradict the ForceCommonNorm above, even though it looks like it should,
    // because CommonNorm applies to the AddPdf while SeparateNorm
    // applies to the ProdPdf.

    comps.clear();
    comps.push_back(bkg3_dalitz);
    // bkg3_dalitz->setDebugMask(1);
    // incsum3->setDebugMask(1);
    comps.push_back(bkg3_dtime);
    // bkg3_dtime->setDebugMask(1);
    // comps.push_back(bkg3_jsugg);
    // sig0_jsugg->setDebugMask(1);

    GooPdf *ret = new ProdPdf("bkg3_total", comps);

    if(fixem)
        ret->setParameterConstantness(true);

    return ret;
}

GooPdf *makeBkg4DalitzPdf(bool fixem = true) {
    if(!kzero_veto)
        makeKzeroVeto();

    comps.clear();
    weights.clear();

    GooPdf *bkg4_dalitz = nullptr;

    if(!notUseBackground4Hist)
        bkg4_dalitz = makeBackgroundHistogram(4);
    else
        bkg4_dalitz = makeBackground4DalitzParam();

    GooPdf *bkg4_dtime = nullptr;

    if(gaussBkgTime)
        bkg4_dtime = makeGaussianTimePdf(4);
    else
        bkg4_dtime = makeExpGausTimePdf(4);

    // Separate sigma_t distribution
    // bkg4_jsugg = makeBkg4_sigma();
    // bkg4_jsugg = sig0_jsugg; // Mikhail uses same sigma distribution for everything.
    bkg4_jsugg = makeBkg_sigma_strips(4);
    bkg4_jsugg->addSpecialMask(PdfBase::ForceSeparateNorm); // See comments to bkg3_jsugg.

    comps.clear();
    comps.push_back(bkg4_dalitz);
    comps.push_back(bkg4_dtime);
    // comps.push_back(bkg4_jsugg);

    ProdPdf *ret = new ProdPdf("bkg4_total", comps);

    if(fixem)
        ret->setParameterConstantness(true);

    return ret;
}

int runCanonicalFit(std::string fname, bool noPlots = true) {
    makeFullFitVariables();

    if(mdslices > 1)
        massd0 = new Observable(
            "massd0", 1.8654 + 0.0075 * md0_lower_window + md0offset, 1.8654 + 0.0075 * md0_upper_window + md0offset);

    std::cout << "Loading events from " << fname << std::endl;
    loadDataFile(fname);

    std::cout << "Creating overall signal PDF\n";
    GooPdf *overallSignal = makeOverallSignal();

    TRandom donram(blindSeed); // The rain and the sun!

    if(0 != blindSeed) {
        dtop0pp._xmixing.setBlind(donram.Gaus(0, 0.005));
        dtop0pp._ymixing.setBlind(donram.Gaus(0, 0.005));
    }

    // overallSignal->setDebugMask(1);

    int oldBins1 = m12->getNumBins();
    int oldBins2 = m13->getNumBins();
    // Too fine a binning here leads to bad results due to fluctuations.
    m12->setNumBins(bkgHistBins);
    m13->setNumBins(bkgHistBins);
    std::cout << "Creating background PDFs\n";
    GooPdf *bkg2Pdf = makeBkg2DalitzPdf();
    GooPdf *bkg3Pdf = makeBkg3DalitzPdf();
    GooPdf *bkg4Pdf = makeBkg4DalitzPdf();
    m12->setNumBins(oldBins1);
    m13->setNumBins(oldBins2);

    getBackgroundFile(2);
    std::cout << "Reading bkg2 parameters from " << strbuffer << std::endl;
    readFromFile(bkg2Pdf, strbuffer);
    // NB, background 3 and 4 params do not actually work. Only hists.
    getBackgroundFile(3);
    std::cout << "Reading bkg3 parameters from " << strbuffer << std::endl;
    readFromFile(bkg3Pdf, strbuffer);
    getBackgroundFile(4);
    std::cout << "Reading bkg4 parameters from " << strbuffer << std::endl;
    readFromFile(bkg4Pdf, strbuffer);

    bkg2Pdf->setParameterConstantness(true);
    bkg3Pdf->setParameterConstantness(true);
    bkg4Pdf->setParameterConstantness(true);

    // bkg3Pdf->setDebugMask(1);

    int eventSize = massd0 ? 11 : 10;
    std::cout << "Setting data size " << eventSize << std::endl;
    signalDalitz->setDataSize(data->getNumEvents(), eventSize); // Must take into account event weights!

    // bkg2Pdf->setDebugMask(1);
    if(incsum1)
        incsum1->setDataSize(data->getNumEvents(), eventSize);

    if(incsum2)
        incsum2->setDataSize(data->getNumEvents(), eventSize);

    // bkg3Pdf->setDebugMask(1);
    if(incsum3)
        incsum3->setDataSize(data->getNumEvents(), eventSize);

    if(incsum4)
        incsum4->setDataSize(data->getNumEvents(), eventSize);

    if(incsum5)
        incsum5->setDataSize(data->getNumEvents(), eventSize);

    if(incsum6)
        incsum6->setDataSize(data->getNumEvents(), eventSize);

    std::cout << "Creating overall PDF\n";
    std::vector<Observable> evtWeights;
    evtWeights.push_back(*wSig0);
    evtWeights.push_back(*wBkg2);
    evtWeights.push_back(*wBkg3);
    evtWeights.push_back(*wBkg4);
    std::vector<PdfBase *> components;
    components.push_back(overallSignal);
    components.push_back(bkg2Pdf);
    components.push_back(bkg3Pdf);
    components.push_back(bkg4Pdf);
    EventWeightedAddPdf *overallPdf = new EventWeightedAddPdf("total", evtWeights, components);
    // overallPdf->setDebugMask(1);
    std::cout << "Copying data to GPU\n";
    overallPdf->setData(data);

    if(paramUp != "") {
        Variable *target = overallPdf->getParameterByName(paramUp);
        assert(target);
        target->setValue(target->getValue() + target->getError());
    }

    if(paramDn != "") {
        Variable *target = overallPdf->getParameterByName(paramDn);
        assert(target);
        target->setValue(target->getValue() - target->getError());
    }

    int retval;
    if(minuit1) {
        GooFit::FitManagerMinuit1 datapdf(overallPdf);
        datapdf.setMaxCalls(64000);
        datapdf.fit();
        retval = datapdf;
    } else {
        GooFit::FitManagerMinuit2 datapdf(overallPdf);
        datapdf.setMaxCalls(64000);
        datapdf.fit();
        retval = datapdf;
    }

#ifdef PROFILING
    overallPdf->printProfileInfo();
#endif

    fmt::print("Fit results Canonical fit:\n"
               "tau    : ({:.3}) fs\n"
               "xmixing: ({:.3})\%\n"
               "ymixing: ({:.3})\%\n",
               1000 * Uncertain(dtop0pp._tau),
               100 * Uncertain(dtop0pp._xmixing),
               100 * Uncertain(dtop0pp._ymixing));

    /*
    std::cout << "Fit results: \n"
        << "tau    : (" << 1000*dtop0pp._tau.getValue() << " $\\pm$ " << 1000*dtop0pp._tau.getError() << ") fs\n"
        << "xmixing: (" << 100*dtop0pp._xmixing.getValue() << " $\\pm$ " << 100*dtop0pp._xmixing.getError() << ")%\n"
        << "ymixing: (" << 100*dtop0pp._ymixing.getValue() << " $\\pm$ " << 100*dtop0pp._ymixing.getError() << ")%\n";
    */

    /*
    double fitx = dtop0pp._xmixing.getValue();
    TH1F xscan("xscan", "", 200, fitx - 0.0001 - 0.0000005, fitx + 0.0001 - 0.0000005);
    xscan.SetStats(false);
    //double fity = dtop0pp._ymixing.getValue();
    for (int i = 0; i < 200; ++i) {
      dtop0pp._xmixing.getValue() = fitx - 0.0001 + 0.000001*i;
      overallPdf->copyParams();
      double nll = overallPdf->calculateNLL();
      printf("%i: %.10f\n", i, nll);
      xscan.SetBinContent(i+1, nll - floor(nll));
    }

    foo->cd();
    xscan.Draw();
    foo->SaveAs("xscan.png");
    */
    if(!noPlots)
        makeDalitzPlots(overallSignal);

    return retval;
}

int runSigmaFit(const char *fname) {
    makeFullFitVariables();

    loM23Sigma = new TH1F("loM23Sigma", "", sigma->getNumBins(), sigma->getLowerLimit(), sigma->getUpperLimit());
    loM23Sigma->SetStats(false);
    hiM23Sigma = new TH1F("hiM23Sigma", "", sigma->getNumBins(), sigma->getLowerLimit(), sigma->getUpperLimit());
    hiM23Sigma->SetStats(false);

    loadDataFile(fname);
    // GooPdf* jsu_gg = makeSignalJSU_gg(-1, false);
    // GooPdf* jsu_gg = makeSigmaMap();
    GooPdf *jsu_gg = makeBkg_sigma_strips(0);
    jsu_gg->setData(data);
    // jsu_gg->copyParams();

    int retval;
    if(minuit1) {
        GooFit::FitManagerMinuit1 datapdf(jsu_gg);
        datapdf.fit();
        retval = datapdf;
    } else {
        GooFit::FitManagerMinuit2 datapdf(jsu_gg);
        datapdf.fit();
        retval = datapdf;
    }

    sprintf(strbuffer, "signal_sigma_%islices_pdf.txt", m23Slices);
    writeToFile(jsu_gg, strbuffer);

    foo->cd();
    plotLoHiSigma();

    std::vector<Observable> gridvars;
    gridvars.push_back(*m12);
    gridvars.push_back(*m13);
    gridvars.push_back(*sigma);
    UnbinnedDataSet grid(gridvars);

    TH1F *sigma_pdfs[6];
    TH1F *sigma_data[6];
    double num_sigma_dat[6];
    double num_sigma_pdf[6];

    for(int i = 0; i < 6; ++i) {
        sprintf(strbuffer, "sigma_pdf_%i", i);
        sigma_pdfs[i] = new TH1F(strbuffer, "", sigma->getNumBins(), sigma->getLowerLimit(), sigma->getUpperLimit());
        sprintf(strbuffer, "sigma_dat_%i", i);
        sigma_data[i] = new TH1F(strbuffer, "", sigma->getNumBins(), sigma->getLowerLimit(), sigma->getUpperLimit());

        num_sigma_dat[i] = 0;
        num_sigma_pdf[i] = 0;

        sigma_data[i]->SetStats(false);
        sigma_data[i]->SetMarkerStyle(8);
        sigma_data[i]->SetMarkerSize(1.2);
        sigma_data[i]->GetXaxis()->SetTitle("Decay time error [ps]");
        sigma_data[i]->GetYaxis()->SetTitle("Events / 8 fs");

        sigma_pdfs[i]->SetStats(false);
        sigma_pdfs[i]->SetLineColor(kBlue);
        sigma_pdfs[i]->SetLineWidth(3);
    }

    double totalPdf = 0;
    double totalDat = 0;

    for(unsigned int evt = 0; evt < data->getNumEvents(); ++evt) {
        double currSigma = data->getValue(*sigma, evt);
        double currm12   = data->getValue(*m12, evt);
        double currm13   = data->getValue(*m13, evt);
        double currm23   = cpuGetM23(currm12, currm13);
        int m23bin       = (int)floor(currm23 / 0.5);
        sigma_data[m23bin]->Fill(currSigma);
        num_sigma_dat[m23bin]++;
        totalDat++;
    }

    for(int i = 0; i < m12->getNumBins(); ++i) {
        m12->setValue(m12->getLowerLimit() + (i + 0.5) * m12->getBinSize());

        for(int j = 0; j < m13->getNumBins(); ++j) {
            m13->setValue(m13->getLowerLimit() + (j + 0.5) * m13->getBinSize());

            if(!cpuDalitz(m12->getValue(), m13->getValue(), _mD0, piZeroMass, piPlusMass, piPlusMass))
                continue;

            for(int k = 0; k < sigma->getNumBins(); ++k) {
                sigma->setValue(sigma->getLowerLimit() + (k + 0.5) * sigma->getBinSize());
                grid.addEvent();
            }
        }
    }

    jsu_gg->setData(&grid);

    std::vector<std::vector<double>> pdfValues = jsu_gg->getCompProbsAtDataPoints();

    for(unsigned int j = 0; j < pdfValues[0].size(); ++j) {
        double currM12   = grid.getValue(*m12, j);
        double currM13   = grid.getValue(*m13, j);
        double currSigma = grid.getValue(*sigma, j);
        double currm23   = cpuGetM23(currM12, currM13);
        int m23bin       = (int)floor(currm23 / 0.5);
        sigma_pdfs[m23bin]->Fill(currSigma, pdfValues[0][j]);
        num_sigma_pdf[m23bin] += pdfValues[0][j];
        totalPdf += pdfValues[0][j];
    }

    for(int i = 1; i <= sigma->getNumBins(); ++i) {
        for(int j = 0; j < 6; ++j) {
            sigma_pdfs[j]->SetBinContent(i, sigma_pdfs[j]->GetBinContent(i) * num_sigma_dat[j] / num_sigma_pdf[j]);
        }
    }

    std::string plotdir = "./plots_from_mixfit/";

    for(int i = 0; i < 6; ++i) {
        if(sigma_data[i]->GetMaximum() > sigma_pdfs[i]->GetMaximum()) {
            sigma_data[i]->Draw("p");
            sigma_pdfs[i]->Draw("lsame");
        } else {
            sigma_pdfs[i]->Draw("l");
            sigma_data[i]->Draw("psame");
        }

        sprintf(strbuffer, "%i", i);
        TText slicenum;
        slicenum.DrawTextNDC(0.2, 0.8, strbuffer);

        foo->SaveAs((plotdir + sigma_pdfs[i]->GetName() + ".png").c_str());
        foo->SetLogy(true);
        foo->SaveAs((plotdir + sigma_pdfs[i]->GetName() + "_log.png").c_str());
        foo->SetLogy(false);
    }

    /*
      // This code assumes you're using the PDF from makeSigmaMap.
    TCanvas foodal("", "", 600, 600);
    foodal.Divide(6, 6, 0, 0);
    for (int i = 0; i < numSigmaBins; ++i) {
      int xbin = i % 6;
      int ybin = i / 6;

      if (0 == sigma_data[i]->numEvents()) continue;

      //m12->getValue() = 0.5*(xbin + 0.5);
      //m13->getValue() = 0.5*(ybin + 0.5);

      std::vector<fptype> values;
      jsuList[i]->evaluateAtPoints(sigma, values);

      double totalPdf = 0;
      for (int bin = 0; bin < sigma->getNumBins(); ++bin) {
        totalPdf += values[bin];
      }
      for (int bin = 0; bin < sigma->getNumBins(); ++bin) {
        sigma_pdf_hists[i]->SetBinContent(bin+1, values[bin] * sigma_data[i]->numEvents() / totalPdf);
      }

      int padNumber = 31 - 6*ybin;
      padNumber += xbin;
      foodal.cd(padNumber);
      sigma_dat_hists[i]->Draw("p");
      sigma_pdf_hists[i]->Draw("lsame");
    }
    foodal.SaveAs("./plots_from_mixfit/sigma_dalitz.png");
    */

    return retval;
}

int runEfficiencyFit(int which) {
    makeFullFitVariables();

    if(3 == which) {
        m12->setNumBins(m12->getNumBins() / 8);
        m13->setNumBins(m13->getNumBins() / 8);
    }

    vector<Observable *> lvars;
    lvars.push_back(m12);
    lvars.push_back(m13);
    // binEffData = new BinnedDataSet(lvars);
    // GooPdf* eff = makeEfficiencyPdf();

    makeKzeroVeto();
    // GooPdf* eff = makeEfficiencyPoly();
    GooPdf *eff = makeEfficiencyPdf();

    std::string fname_3flat = app_ptr->get_filename("dataFiles/efficiency_bkg3_flat.txt", "examples/pipipi0DPFit");
    std::string fname_flat  = app_ptr->get_filename("dataFiles/efficiency_flat.txt", "examples/pipipi0DPFit");

    if(3 == which)
        loadDataFile(fname_3flat);
    else
        loadDataFile(fname_flat);

    if(underlyingBins) {
        underlyingBins->GetZaxis()->SetRangeUser(10, 40);
        underlyingBins->Draw("colz");
        foo->SaveAs("./plots_from_mixfit/efficiency_bins.png");
    }

    // eff->setDebugMask(1);
    eff->setData(data);

    int retval;
    if(minuit1) {
        GooFit::FitManagerMinuit1 datapdf(eff);
        datapdf.fit();
        retval = datapdf;
    } else {
        GooFit::FitManagerMinuit2 datapdf(eff);
        datapdf.fit();
        retval = datapdf;
    }

    // plotFit(sigma, data, jsu_gg);

    TH2F dalitz_dat_hist("dalitz_dat_hist",
                         "",
                         m12->getNumBins(),
                         m12->getLowerLimit(),
                         m12->getUpperLimit(),
                         m13->getNumBins(),
                         m13->getLowerLimit(),
                         m13->getUpperLimit());
    dalitz_dat_hist.SetStats(false);
    dalitz_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{0}) [GeV]");
    dalitz_dat_hist.GetYaxis()->SetTitle("m^{2}(#pi^{-} #pi^{0}) [GeV]");
    TH2F dalitz_pdf_hist("dalitz_pdf_hist",
                         "",
                         m12->getNumBins(),
                         m12->getLowerLimit(),
                         m12->getUpperLimit(),
                         m13->getNumBins(),
                         m13->getLowerLimit(),
                         m13->getUpperLimit());
    dalitz_pdf_hist.SetStats(false);

    double totalPdf = 0;
    double totalDat = 0;

    for(unsigned int evt = 0; evt < data->getNumEvents(); ++evt) {
        double currval = data->getValue(*m12, evt);
        // m12_dat_hist.Fill(currval);
        double currval2 = data->getValue(*m13, evt);
        // m13_dat_hist.Fill(currval2);
        dalitz_dat_hist.Fill(currval, currval2);
        totalDat++;
    }

    std::vector<Observable> nvars;
    nvars.push_back(*m12);
    nvars.push_back(*m13);
    UnbinnedDataSet currData(nvars);

    for(int i = 0; i < m12->getNumBins(); ++i) {
        m12->setValue(m12->getLowerLimit() + (i + 0.5) * m12->getBinSize());

        for(int j = 0; j < m13->getNumBins(); ++j) {
            m13->setValue(m13->getLowerLimit() + (j + 0.5) * m13->getBinSize());

            if(!cpuDalitz(m12->getValue(), m13->getValue(), _mD0, piZeroMass, piPlusMass, piPlusMass))
                continue;

            currData.addEvent();
        }
    }

    eff->setData(&currData);

    // eff->setDebugMask(1);
    std::vector<std::vector<double>> pdfValues = eff->getCompProbsAtDataPoints();

    for(unsigned int j = 0; j < pdfValues[0].size(); ++j) {
        double currVal = currData.getValue(*m12, j);
        // m12_pdf_hist.Fill(currVal, pdfValues[0][j]);
        double currVal2 = currData.getValue(*m13, j);
        // m13_pdf_hist.Fill(currVal, pdfValues[0][j]);
        dalitz_pdf_hist.Fill(currVal, currVal2, pdfValues[0][j]);

        totalPdf += pdfValues[0][j];

        if(std::isnan(pdfValues[0][j])) {
            std::cout << "Major problem: " << currVal << " " << currVal2 << " " << j << std::endl;
            assert(false);
        }

        if(std::isinf(pdfValues[0][j])) {
            std::cout << "Infinity " << j << std::endl;
            assert(false);
        }
    }

    /*for (int i = 1; i <= m12->getNumBins(); ++i) {
      m12_pdf_hist.SetBinContent(i, m12_pdf_hist.GetBinContent(i) * totalDat / totalPdf);
    }

    for (int i = 1; i <= m13->getNumBins(); ++i) {
      m13_pdf_hist.SetBinContent(i, m13_pdf_hist.GetBinContent(i) * totalDat / totalPdf);
    }
    */
    for(int i = 1; i <= m12->getNumBins(); ++i) {
        for(int j = 1; j <= m13->getNumBins(); ++j) {
            double currNormEff = dalitz_pdf_hist.GetBinContent(i, j) * totalDat / totalPdf;
            dalitz_pdf_hist.SetBinContent(i, j, currNormEff);

            if(currNormEff >= 50)
                std::cout << "High efficiency: " << currNormEff << " " << i << " " << j << std::endl;
        }
    }

    foodal->cd();

    dalitz_dat_hist.GetZaxis()->SetRangeUser(0, 40);
    dalitz_dat_hist.Draw("colz");
    foodal->SaveAs("./plots_from_mixfit/efficiency_dat.png");

    dalitz_pdf_hist.GetZaxis()->SetRangeUser(0, 40);
    dalitz_pdf_hist.Draw("colz");
    foodal->SaveAs("./plots_from_mixfit/efficiency_pdf.png");
    foo->cd();

    TH1F pullplot("pullplot", "", 100, -5, 5);
    TH1F hiM23pullplot("hiM23pullplot", "", 100, -5, 5);
    TH1F loM23pullplot("loM23pullplot", "", 100, -5, 5);

    for(int i = 1; i <= m12->getNumBins(); ++i) {
        double m12loedge
            = m12->getLowerLimit() + ((m12->getUpperLimit() - m12->getLowerLimit()) / m12->getNumBins()) * (i - 1);
        double m12hiedge
            = m12->getLowerLimit() + ((m12->getUpperLimit() - m12->getLowerLimit()) / m12->getNumBins()) * (i);

        for(int j = 1; j <= m13->getNumBins(); ++j) {
            double m13loedge
                = m13->getLowerLimit() + ((m13->getUpperLimit() - m13->getLowerLimit()) / m13->getNumBins()) * (j - 1);

            if(!cpuDalitz(m12loedge, m13loedge, _mD0, piZeroMass, piPlusMass, piPlusMass)) {
                dalitz_dat_hist.SetBinContent(i, j, 0);
                continue;
            }

            if(!cpuDalitz(m12hiedge, m13loedge, _mD0, piZeroMass, piPlusMass, piPlusMass)) {
                dalitz_dat_hist.SetBinContent(i, j, 0);
                continue;
            }

            double m13hiedge
                = m13->getLowerLimit() + ((m13->getUpperLimit() - m13->getLowerLimit()) / m13->getNumBins()) * (j);

            if(!cpuDalitz(m12loedge, m13hiedge, _mD0, piZeroMass, piPlusMass, piPlusMass)) {
                dalitz_dat_hist.SetBinContent(i, j, 0);
                continue;
            }

            if(!cpuDalitz(m12hiedge, m13hiedge, _mD0, piZeroMass, piPlusMass, piPlusMass)) {
                dalitz_dat_hist.SetBinContent(i, j, 0);
                continue;
            }

            double dat = dalitz_dat_hist.GetBinContent(i, j);
            double pdf = dalitz_pdf_hist.GetBinContent(i, j);

            double pull = (dat - pdf) / sqrt(max(1.0, dat));
            // if (fabs(pull) > 5) continue;
            dalitz_dat_hist.SetBinContent(i, j, pull);
            pullplot.Fill(pull);

            double currm12 = dalitz_dat_hist.GetXaxis()->GetBinCenter(i);
            double currm13 = dalitz_dat_hist.GetYaxis()->GetBinCenter(j);

            // double currm23 = cpuGetM23(currm12, currm13);
            // if (currm23 > 1.5) hiM23pullplot.Fill(pull);
            // else loM23pullplot.Fill(pull);
            if((currm13 > 2) || (currm12 > 2))
                hiM23pullplot.Fill(pull);
            else
                loM23pullplot.Fill(pull);
        }
    }

    foodal->cd();
    dalitz_dat_hist.GetZaxis()->SetRangeUser(-5, 5);
    dalitz_dat_hist.Draw("colz");
    foodal->SaveAs("./plots_from_mixfit/efficiency_pull.png");

    foo->cd();
    pullplot.Draw();
    foo->SaveAs("./plots_from_mixfit/effpull.png");

    hiM23pullplot.Draw();
    foo->SaveAs("./plots_from_mixfit/hieffpull.png");

    loM23pullplot.Draw();
    foo->SaveAs("./plots_from_mixfit/loeffpull.png");

    return retval;
}

int runBackgroundDalitzFit(int bkgType, bool plots) {
    makeFullFitVariables();
    makeKzeroVeto();

    GooPdf *bkgPdf = 0;

    switch(bkgType) {
    default:
    case 2:
        bkgPdf = makeBkg2DalitzPdf(false);
        break;

    case 3:
        bkgPdf = makeBkg3DalitzPdf(false);
        break;

    case 4:
        bkgPdf = makeBkg4DalitzPdf(false);
        break;
    }

    sprintf(strbuffer, "./dataFiles/bkgDalitz_%i.txt", bkgType);
    std::string fname = app_ptr->get_filename(strbuffer, "examples/pipipi0DPFit");
    loadDataFile(fname);

    bkgPdf->setData(data);
    // bkgPdf->setDebugMask(1);

    // Incoherent-sum components need to know data size; check which ones exist.
    int eventSize = 5; // m12 m13 dtime sigma eventNumber

    if(incsum1)
        incsum1->setDataSize(data->getNumEvents(), eventSize);

    if(incsum2)
        incsum2->setDataSize(data->getNumEvents(), eventSize);

    if(incsum3)
        incsum3->setDataSize(data->getNumEvents(), eventSize);

    if(incsum4)
        incsum4->setDataSize(data->getNumEvents(), eventSize);

    if(incsum5)
        incsum5->setDataSize(data->getNumEvents(), eventSize);

    if(incsum6)
        incsum6->setDataSize(data->getNumEvents(), eventSize);

    int retval;
    if(minuit1) {
        GooFit::FitManagerMinuit1 fitter(bkgPdf);
        fitter.setMaxCalls(32000);
        fitter.fit();
        retval = fitter;
    } else {
        GooFit::FitManagerMinuit2 fitter(bkgPdf);
        fitter.setMaxCalls(32000);
        fitter.fit();
        retval = fitter;
    }

    if(plots) {
        // sprintf(strbuffer, "./plots_from_mixfit/bkgdalitz_%i/", bkgType);
        // makeDalitzPlots(bkgPdf, strbuffer);
        getBackgroundFile(bkgType);
        std::string fname = app_ptr->get_filename(strbuffer, "examples/pipipi0DPFit");
        writeToFile(bkgPdf, fname.c_str());
    }

    return retval;
}

void getBackgroundFile(int bkgType) {
    if(mikhailSetup)
        sprintf(strbuffer, "./bkg_%i_mikhail.txt", bkgType);
    else {
        if(2 == bkgType) {
            if(Sideband == bkg2Model)
                sprintf(strbuffer, "./bkg_2_pdf_sideband_%islices.txt", m23Slices);
            else
                sprintf(strbuffer, "./bkg_2_pdf_%islices.txt", m23Slices);
        } else {
            std::string pdftype;

            if(((3 == bkgType) && (notUseBackground3Hist)) || ((4 == bkgType) && (notUseBackground4Hist)))
                pdftype = "_param";

            sprintf(strbuffer, "./bkg_%i_pdf%s.txt", bkgType, pdftype.c_str());
        }
    }
}

void makeTimePlots(std::string fname) {
    makeFullFitVariables();
    massd0 = new Observable("massd0", 1.8654 + 0.0075 * md0_lower_window, 1.8654 + 0.0075 * md0_upper_window);
    massd0->setNumBins(180);
    std::cout << "Loading MC data from " << fname << std::endl;
    loadDataFile(fname);

    TH1F timeMean("timeMean", "", 6, massd0->getLowerLimit(), massd0->getUpperLimit());
    timeMean.SetStats(false);
    timeMean.SetLineWidth(3);
    timeMean.SetXTitle("#pi#pi#pi^{0} mass [GeV]");
    timeMean.SetYTitle("Mean of decay time [ps]");
    TH2F timeVsMass("timeVsMass",
                    "",
                    massd0->getNumBins(),
                    massd0->getLowerLimit(),
                    massd0->getUpperLimit(),
                    dtime->getNumBins(),
                    dtime->getLowerLimit(),
                    dtime->getUpperLimit());
    timeVsMass.SetStats(false);
    timeVsMass.GetXaxis()->SetTitle("#pi#pi#pi^{0} mass [GeV]");
    timeVsMass.GetYaxis()->SetTitle("Decay time [ps]");

    int colors[6] = {kViolet + 1, kBlue, kCyan, kGreen, kYellow, kRed};
    TH1F *timePlots[6];
    TH1F *massPlots[5];

    for(int i = 0; i < 6; ++i) {
        sprintf(strbuffer, "timePlot_%i.png", i);
        timePlots[i] = new TH1F(strbuffer, "", dtime->getNumBins(), dtime->getLowerLimit(), dtime->getUpperLimit());
        timePlots[i]->SetStats(false);
        timePlots[i]->SetXTitle("Decay time [ps]");
        timePlots[i]->SetYTitle("Ratio");
        timePlots[i]->SetLineWidth(3);
        timePlots[i]->SetLineColor(colors[i]);

        if(i == 5)
            continue;

        sprintf(strbuffer, "massPlot_%i.png", i);
        massPlots[i] = new TH1F(strbuffer, "", massd0->getNumBins(), massd0->getLowerLimit(), massd0->getUpperLimit());
        massPlots[i]->SetStats(false);
        massPlots[i]->SetLineWidth(3);
        massPlots[i]->SetLineColor(colors[i]);
    }

    for(int i = 0; i < data->getNumEvents(); ++i) {
        data->loadEvent(i);
        timeVsMass.Fill(massd0->getValue(), dtime->getValue());

        if(massd0->getValue() >= massd0->getUpperLimit())
            continue;

        if(massd0->getValue() < massd0->getLowerLimit())
            continue;

        int slice = (int)floor(6 * (massd0->getValue() - massd0->getLowerLimit())
                               / (massd0->getUpperLimit() - massd0->getLowerLimit()));
        timePlots[slice]->Fill(dtime->getValue());

        slice = (int)floor(5 * (dtime->getValue() - dtime->getLowerLimit())
                           / (dtime->getUpperLimit() - dtime->getLowerLimit()));
        massPlots[slice]->Fill(massd0->getValue());
    }

    foo->cd();
    normalize(timePlots[3]);
    timePlots[3]->SetMinimum(0);
    timePlots[3]->Draw("hist");

    for(int i = 0; i < 6; ++i) {
        normalize(timePlots[i]);
        timePlots[i]->Draw("histsame");
        timeMean.SetBinContent(i + 1, timePlots[i]->GetMean());
        timeMean.SetBinError(i + 1, timePlots[i]->GetMeanError());
    }

    foo->SaveAs("timePlots.png");
    timeMean.Draw("e");
    foo->SaveAs("timeMeanPlot.png");

    // normalize(massPlots[2]);
    massPlots[2]->GetYaxis()->SetRangeUser(0, massPlots[2]->GetMaximum() * 1.1);
    massPlots[2]->Draw("");

    for(int i = 0; i < 5; ++i) {
        // normalize(massPlots[i]);
        massPlots[i]->Draw("same");
    }

    foo->SaveAs("massPlots.png");

    timeVsMass.Draw("colz");

    for(int i = 0; i < 6; ++i) {
        TLine *currLine
            = new TLine(massd0->getLowerLimit() + (i + 0) * (massd0->getUpperLimit() - massd0->getLowerLimit()) / 6,
                        dtime->getLowerLimit() + 0.09,
                        massd0->getLowerLimit() + (i + 1) * (massd0->getUpperLimit() - massd0->getLowerLimit()) / 6,
                        dtime->getLowerLimit() + 0.09);
        currLine->SetLineWidth(12);
        currLine->SetLineColor(colors[i]);
        currLine->Draw();

        if(5 == i)
            continue;

        currLine = new TLine(massd0->getLowerLimit() + 0.00025,
                             dtime->getLowerLimit() + (i + 0) * (dtime->getUpperLimit() - dtime->getLowerLimit()) / 5,
                             massd0->getLowerLimit() + 0.00025,
                             dtime->getLowerLimit() + (i + 1) * (dtime->getUpperLimit() - dtime->getLowerLimit()) / 5);
        currLine->SetLineWidth(12);
        currLine->SetLineColor(colors[i]);
        //    currLine->Draw();
    }

    foo->SaveAs("timeVsMass.png");
}

int runBackgroundSigmaFit(int bkgType) {
    makeFullFitVariables();

    GooPdf *bkgPdf = 0;

    switch(bkgType) {
    default:
    case 2:
        bkgPdf = makeBkg2_sigma();
        break;

    case 3:
        loM23Sigma = new TH1F("loM23Sigma", "", sigma->getNumBins(), sigma->getLowerLimit(), sigma->getUpperLimit());
        loM23Sigma->SetStats(false);
        hiM23Sigma = new TH1F("hiM23Sigma", "", sigma->getNumBins(), sigma->getLowerLimit(), sigma->getUpperLimit());
        hiM23Sigma->SetStats(false);

        bkgPdf = makeBkg_sigma_strips(3);
        break;
    }

    sprintf(strbuffer, "./dataFiles/bkgDalitz_%i.txt", bkgType);
    std::string fname = app_ptr->get_filename(strbuffer, "examples/pipipi0DPFit");
    loadDataFile(fname);
    bkgPdf->setData(data);

    int retval;

    if(minuit1) {
        GooFit::FitManagerMinuit1 fitter(bkgPdf);
        fitter.setMaxCalls(8000);
        fitter.fit();
        retval = fitter;
    } else {
        GooFit::FitManagerMinuit2 fitter(bkgPdf);
        fitter.setMaxCalls(8000);
        fitter.fit();
        retval = fitter;
    }

    // bkgPdf->setDebugMask(1);
    plotFit(*sigma, data, bkgPdf);
    plotLoHiSigma();

    return retval;

    // sprintf(strbuffer, "./plots_from_mixfit/bkgdalitz_%i/", bkgType);
    // makeDalitzPlots(bkgPdf, strbuffer);
}

void writeBackgroundHistograms(int bkg) {
    makeFullFitVariables();
    SmoothHistogramPdf *bkg3 = makeBackgroundHistogram(bkg);
    thrust::host_vector<fptype> host_hist;
    bkg3->extractHistogram(host_hist);
    sprintf(strbuffer, "Bkg%i_dalitzhist.txt", bkg);
    std::string fname = app_ptr->get_filename(strbuffer, "examples/pipipi0DPFit");
    writeListOfNumbers(host_hist, fname.c_str());
}

void set_bkg_model_from_string() {
    if(bkg2Model_str == "histogram")
        bkg2Model = Histogram;
    else if(bkg2Model_str == "parameter")
        bkg2Model = Parameter;
    else if(bkg2Model_str == "sideband")
        bkg2Model = Sideband;

    if(mikhailSetup) {
        m23Slices    = 1;
        gaussBkgTime = true;
    }
}

void parseArg(GooFit::App *app) {
    app->add_option("--luckyFrac", luckyFrac, "", true);
    app->add_option("--mesonRad", mesonRad, "", true);
    app->add_option("--normBins", normBinning, "", true);
    app->add_option("--blindSeed", blindSeed, "", true);
    app->add_option("--mdslices", mdslices, "", true);
    app->add_option("--offset", md0offset, "Offest in GeV", true);
    // Previously in MeV
    app->add_option("--upper_window", md0_upper_window, "", true);
    app->add_option("--lower_window", md0_lower_window, "", true);
    app->add_option("--upper_delta_window", deltam_upper_window, "", true);
    app->add_option("--lower_delta_window", deltam_lower_window, "", true);
    app->add_option("--upperTime", upperTime, "", true);
    app->add_option("--lowerTime", lowerTime, "", true);
    app->add_option("--maxSigma", maxSigma, "", true);
    app->add_option("--polyEff", polyEff, "", true);
    app->add_option("--m23Slices", m23Slices, "", true);
    app->add_option("--bkgRandSeed", bkgHistRandSeed, "", true);

    app->add_flag("--drop-rho_1450", drop_rho_1450);
    app->add_flag("--drop-rho_1700", drop_rho_1700);
    app->add_flag("--drop-f0_980", drop_f0_980);
    app->add_flag("--drop-f0_1370", drop_f0_1370);
    app->add_flag("--drop-f0_1500", drop_f0_1500);
    app->add_flag("--drop-f0_1710", drop_f0_1710);
    app->add_flag("--drop-f2_1270", drop_f2_1270);
    app->add_flag("--drop-f0_600", drop_f0_600);

    app->add_flag("--histSigma", useHistogramSigma);
    app->add_flag("--makePlots", makePlots);
    app->add_set("--mkg2Model", bkg2Model_str, {"histogram", "parameter", "sideband"}, "", true);
    app->add_flag("--bkg3Hist", notUseBackground3Hist);
    app->add_flag("--bkg4Hist", notUseBackground4Hist);
    app->add_option("--bkgHistBins", bkgHistBins, "", true);
    app->add_option("--varyParameterUp", paramUp, "", true);
    app->add_option("--varyParameterDn", paramDn, "", true);
    app->add_flag("--mikhail", mikhailSetup);
}

int main(int argc, char **argv) {
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

    int retval = 0;

    GooFit::Application app("pipipi0 Dalitz fit example", argc, argv);
    app_ptr = &app;
    app.require_subcommand();

    app.add_flag("--minuit1", minuit1, "Use Minuit 1 instead of Minuit 2");

    std::string data;
    int sample = 0;
    int load   = 1;
    bool plots;
    int genResolutions = 0;
    double dplotres    = 0;

    auto toy = app.add_subcommand("toy", "Toy MC Performance evaluation");
    toy->add_option("-s,--sample,sample", sample, "Sample number to use", true);
    toy->add_option("-l,--load,load", load, "Number of times to load", true);
    toy->add_option("-m,--max", maxEvents, "Maximum number of events to read", true);
    toy->add_flag("-p,--plot", plots, "Also make plots");
    toy->set_callback([&]() { retval = runToyFit(sample, load, plots); });

    auto truth_fit = app.add_subcommand("truth", "Truth Monte Carlo fit");
    truth_fit->add_option("-d,--data,data", data, "Data to use")->required()->check(GooFit::ExistingFile);
    truth_fit->set_callback([&]() { retval = runTruthMCFit(data, false); });

    auto sigma_fit = app.add_subcommand("sigma", "Run sigma fit");
    sigma_fit->add_option("-d,--data,data", data, "Data to use")->required()->check(GooFit::ExistingFile);
    sigma_fit->add_option("-s,--slices,slices", m23Slices, "m23 slices")->required();
    sigma_fit->set_callback([&]() { retval = runSigmaFit(data.c_str()); });

    auto efficiency_fit = app.add_subcommand("efficiency", "Run efficiency fit");
    efficiency_fit->add_option("-s,--sample,sample", sample, "Sample number to use", true);
    efficiency_fit->set_callback([&]() { retval = runEfficiencyFit(sample); });

    auto canonical_fit = app.add_subcommand("canonical", "Run the canonical fit");
    canonical_fit->add_option("-d,--data,data", data, "Data to use")->required()->check(GooFit::ExistingFile);
    parseArg(canonical_fit);
    canonical_fit->set_callback([&]() {
        set_bkg_model_from_string();
        retval = runCanonicalFit(data, !makePlots);
    });

    auto background_dalitz_fit = app.add_subcommand("background_dalitz", "Run the background Dalitz fit");
    background_dalitz_fit->add_option("-s,--sample,sample", sample, "Sample number to use", true);
    parseArg(background_dalitz_fit);
    background_dalitz_fit->set_callback([&]() {
        set_bkg_model_from_string();
        retval = runBackgroundDalitzFit(sample, true);
    });

    auto background_sigma_fit = app.add_subcommand("background_sigma", "Run background sigma fit");
    background_sigma_fit->add_option("-s,--sample,sample", sample, "Sample number to use", true);
    background_sigma_fit->set_callback([&]() { retval = runBackgroundSigmaFit(sample); });

    auto write_background_histograms = app.add_subcommand("background_histograms", "Write background histograms");
    write_background_histograms->add_option("-s,--sample,sample", sample, "Sample number to use", true);
    write_background_histograms->set_callback([&]() { writeBackgroundHistograms(sample); });

    auto run_gen_mc_fit = app.add_subcommand("run_gen_mc", "Run generated Monte Carlo fit");
    run_gen_mc_fit->add_option("-d,--data,data", data, "Data to use")->required()->check(GooFit::ExistingFile);
    run_gen_mc_fit->add_option("-g,--genres,gen-resolutions", genResolutions)->required();
    run_gen_mc_fit->add_option("-p,--dplotres,dplotres", dplotres);
    run_gen_mc_fit->set_callback([&]() {
        if(!(DplotRes & genResolutions))
            dplotres = 0;
        retval = runGeneratedMCFit(data, genResolutions, dplotres);
    });

    auto make_time_plots = app.add_subcommand("make_time_plots", "Make time plots");
    make_time_plots->add_option("-d,--data,data", data, "Data to use")->required()->check(GooFit::ExistingFile);
    make_time_plots->set_callback([&]() { makeTimePlots(data); });

    GOOFIT_PARSE(app);

    return retval;
}
