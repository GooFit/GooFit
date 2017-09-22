#include "goofit/Application.h"
#include "goofit/Variable.h"
#include "goofit/FitManager.h"
#include "goofit/BinnedDataSet.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/FitControl.h"

#include "goofit/PDFs/basic/PolynomialPdf.h"
#include "TMinuit.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLatex.h"

#include <sys/time.h>
#include <sys/times.h>

TCanvas foo;
timeval startTime, stopTime, totalTime;

#include <vector>
#include <iostream>
#include <string>

using namespace std;
using namespace GooFit;

Variable *decayTime  = 0;
Variable *constaCoef = 0;
Variable *linearCoef = 0;
Variable *secondCoef = 0;

double integralExpCon(double lo, double hi) { return (exp(-lo) - exp(-hi)); }

double integralExpLin(double lo, double hi) { return ((lo + 1) * exp(-lo) - (hi + 1) * exp(-hi)); }

double integralExpSqu(double lo, double hi) {
    return ((lo * lo + 2 * lo + 2) * exp(-lo) - (hi * hi + 2 * hi + 2) * exp(-hi));
}

void generateEvents(vector<int> &rsEvtVec,
                    vector<int> &wsEvtVec,
                    Variable const *const decayTime,
                    double conCoef,
                    double linCoef,
                    double squCoef,
                    int eventsToGenerate) {
    static TRandom donram(24);
    double totalRSintegral = integralExpCon(0, 100);
    double step            = (decayTime->getUpperLimit() - decayTime->getLowerLimit()) / decayTime->getNumBins();

    for(int i = 0; i < decayTime->getNumBins(); ++i) {
        double binStart = i * step;
        binStart += decayTime->getLowerLimit();
        double binFinal = binStart + step;

        double rsIntegral = integralExpCon(binStart, binFinal);
        double wsIntegral = conCoef * integralExpCon(binStart, binFinal);
        wsIntegral += linCoef * integralExpLin(binStart, binFinal);
        wsIntegral += squCoef * integralExpSqu(binStart, binFinal);

        double expectedRSevts = eventsToGenerate * rsIntegral / totalRSintegral;
        double expectedWSevts = eventsToGenerate * wsIntegral / totalRSintegral;

        int rsEvts  = donram.Poisson(expectedRSevts);
        int wsEvts  = donram.Poisson(expectedWSevts);
        rsEvtVec[i] = rsEvts;
        wsEvtVec[i] = wsEvts;

        if(0 == (i % 10))
            std::cout << "Events in bin " << i << " : " << rsEvts << " (" << expectedRSevts << ") " << wsEvts << " ("
                      << expectedWSevts << ")\n";
    }
}

int fitRatio(vector<int> &rsEvts, vector<int> &wsEvts, std::string plotName = "") {
    TH1D *ratioHist
        = new TH1D("ratioHist", "", decayTime->getNumBins(), decayTime->getLowerLimit(), decayTime->getUpperLimit());

    BinnedDataSet *ratioData = new BinnedDataSet(decayTime);

    for(unsigned int i = 0; i < wsEvts.size(); ++i) {
        double ratio = wsEvts[i];

        if(0 == rsEvts[i])
            rsEvts[i] = 1; // Cheating to avoid div by zero.

        ratio /= rsEvts[i];

        if(0 == wsEvts[i])
            wsEvts[i] = 1; // Avoid zero errors

        double error = wsEvts[i] / pow(rsEvts[i], 2);
        error += pow(wsEvts[i], 2) / pow(rsEvts[i], 3);
        error = sqrt(error);

        ratioData->setBinContent(i, ratio);
        ratioData->setBinError(i, error);
        ratioHist->SetBinContent(i + 1, ratio);
        ratioHist->SetBinError(i + 1, error);
    }

    if(0 == constaCoef) {
        constaCoef = new Variable("constaCoef", 0.03, 0.01, -1, 1);
        constaCoef->setValue(0.03);
        constaCoef->setError(0.01);
        linearCoef = new Variable("linearCoef", 0, 0.01, -1, 1);
        linearCoef->setValue(0.00);
        linearCoef->setError(0.01);
        secondCoef = new Variable("secondCoef", 0, 0.01, -1, 1);
        secondCoef->setValue(0.00);
        secondCoef->setError(0.01);
    }

    vector<Variable *> weights;
    weights.push_back(constaCoef);
    weights.push_back(linearCoef);
    weights.push_back(secondCoef);

    PolynomialPdf *poly = new PolynomialPdf("poly", decayTime, weights);
    poly->setFitControl(new BinnedErrorFit());
    poly->setData(ratioData);
    FitManager datapdf{poly};

    gettimeofday(&startTime, nullptr);
    datapdf.fit();
    gettimeofday(&stopTime, nullptr);

    vector<fptype> values = poly->evaluateAtPoints(decayTime);
    TH1D pdfHist("pdfHist", "", decayTime->getNumBins(), decayTime->getLowerLimit(), decayTime->getUpperLimit());

    for(int i = 0; i < values.size(); ++i) {
        pdfHist.SetBinContent(i + 1, values[i]);
    }

    ratioHist->SetMarkerStyle(8);
    ratioHist->SetMarkerSize(0.5);
    ratioHist->SetStats(false);
    ratioHist->Draw("p");

    char strbuffer[1000];
    sprintf(
        strbuffer, "Constant [10^{-2}] : %.3f #pm %.3f", 1e2 * constaCoef->getValue(), constaCoef->getError() * 1e2);
    TLatex res1(0.14, 0.83, strbuffer);
    res1.SetNDC(true);
    sprintf(
        strbuffer, "Linear [10^{-4}]   : %.3f #pm %.3f", 1e4 * linearCoef->getValue(), linearCoef->getError() * 1e4);
    TLatex res2(0.14, 0.73, strbuffer);
    res2.SetNDC(true);
    sprintf(
        strbuffer, "Quadratic [10^{-6}]: %.3f #pm %.3f", 1e6 * secondCoef->getValue(), secondCoef->getError() * 1e6);
    TLatex res3(0.14, 0.63, strbuffer);
    res3.SetNDC(true);

    res1.Draw();
    res2.Draw();
    res3.Draw();

    pdfHist.SetLineColor(kBlue);
    pdfHist.SetLineWidth(3);
    pdfHist.SetStats(false);
    pdfHist.Draw("lsame");
    foo.SaveAs(plotName.c_str());

    std::cout << "Polynomial function: " << poly->getCoefficient(2) << " * t^2 + " << poly->getCoefficient(1)
              << " * t + " << poly->getCoefficient(0) << std::endl;

    delete ratioHist;
    delete ratioData;
    delete poly;

    return datapdf;
}

double dzero_con     = 0;
double dzero_lin     = 0;
double dzero_qua     = 0;
double dzero_con_err = 0;
double dzero_lin_err = 0;
double dzero_qua_err = 0;
double d0bar_con     = 0;
double d0bar_lin     = 0;
double d0bar_qua     = 0;
double d0bar_con_err = 0;
double d0bar_lin_err = 0;
double d0bar_qua_err = 0;

vector<double> ratios;
vector<double> errors;

void cpvFitFcn(int &npar, double *gin, double &fun, double *fp, int iflag) {
    double conCoef = fp[0];
    double linCoef = fp[1];
    double squCoef = fp[2];

    double chisq = 0;
    double step  = (decayTime->getUpperLimit() - decayTime->getLowerLimit()) / decayTime->getNumBins();

    for(unsigned int i = 0; i < ratios.size(); ++i) {
        double currDTime = decayTime->getLowerLimit() + (i + 0.5) * step;
        double pdfval    = conCoef + linCoef * currDTime + squCoef * currDTime * currDTime;
        chisq += pow((pdfval - ratios[i]) / errors[i], 2);
    }

    fun = chisq;
}

void fitRatioCPU(vector<int> &rsEvts, vector<int> &wsEvts) {
    TH1D *ratioHist
        = new TH1D("ratioHist", "", decayTime->getNumBins(), decayTime->getLowerLimit(), decayTime->getUpperLimit());

    ratios.resize(wsEvts.size());
    errors.resize(wsEvts.size());

    for(unsigned int i = 0; i < wsEvts.size(); ++i) {
        double ratio = wsEvts[i];

        if(0 == rsEvts[i])
            rsEvts[i] = 1; // Cheating to avoid div by zero.

        ratio /= rsEvts[i];

        if(0 == wsEvts[i])
            wsEvts[i] = 1; // Avoid zero errors

        double error = wsEvts[i] / pow(rsEvts[i], 2);
        error += pow(wsEvts[i], 2) / pow(rsEvts[i], 3);
        error = sqrt(error);

        ratios[i] = ratio;
        errors[i] = error;
        ratioHist->SetBinContent(i + 1, ratio);
        ratioHist->SetBinError(i + 1, error);
    }

    TMinuit *minuit = new TMinuit(3);
    minuit->DefineParameter(0, "constaCoef", 0.03, 0.01, -1, 1);
    minuit->DefineParameter(1, "linearCoef", 0, 0.01, -1, 1);
    minuit->DefineParameter(2, "secondCoef", 0, 0.01, -1, 1);
    minuit->SetFCN(cpvFitFcn);

    gettimeofday(&startTime, nullptr);
    minuit->Migrad();
    gettimeofday(&stopTime, nullptr);
}

int main(int argc, char **argv) {
    GooFit::Application app("Chi-square example", argc, argv);

    int numbins = 100;
    app.add_option("-n,--numbins", numbins, "Number of bins", true);

    try {
        app.run();
    } catch(const GooFit::ParseError &e) {
        return app.exit(e);
    }

    // Time is in units of lifetime
    Variable decayTime{"decayTime", 100, 0, 10};
    decayTime.setNumBins(numbins);

    double rSubD = 0.03;
    double rBarD = 0.03;
    double delta = 0;
    double wpPhi = 0;
    double x_mix = 0.0016;
    double y_mix = 0.0055;
    double magPQ = 1.0;
    double magQP = 1.0 / magPQ;

    int eventsToGenerate = 10000000;

    vector<int> dZeroEvtsWS(decayTime.getNumBins());
    vector<int> dZeroEvtsRS(decayTime.getNumBins());
    vector<int> d0barEvtsWS(decayTime.getNumBins());
    vector<int> d0barEvtsRS(decayTime.getNumBins());

    double dZeroLinearCoef = magPQ * sqrt(rSubD) * (y_mix * cos(delta + wpPhi) - x_mix * sin(delta + wpPhi));
    double d0barLinearCoef = magQP * sqrt(rBarD) * (y_mix * cos(delta - wpPhi) - x_mix * sin(delta - wpPhi));

    double dZeroSecondCoef = 0.25 * magPQ * magPQ * (x_mix * x_mix + y_mix * y_mix);
    double d0barSecondCoef = 0.25 * magQP * magQP * (x_mix * x_mix + y_mix * y_mix);

    generateEvents(dZeroEvtsRS, dZeroEvtsWS, &decayTime, rSubD, dZeroLinearCoef, dZeroSecondCoef, eventsToGenerate);
    generateEvents(d0barEvtsRS, d0barEvtsWS, &decayTime, rBarD, d0barLinearCoef, d0barSecondCoef, eventsToGenerate);

    double gpuTime = 0;
    double cpuTime = 0;

    int retval;
    retval = fitRatio(dZeroEvtsRS, dZeroEvtsWS, "dzeroEvtRatio.png");

    if(retval != 0)
        return retval;

    timersub(&stopTime, &startTime, &totalTime);
    gpuTime += totalTime.tv_sec + totalTime.tv_usec / 1000000.0;
    retval = fitRatio(d0barEvtsRS, d0barEvtsWS, "dzbarEvtRatio.png");
    if(retval != 0)
        return retval;
    timersub(&stopTime, &startTime, &totalTime);
    gpuTime += totalTime.tv_sec + totalTime.tv_usec / 1000000.0;

    fitRatioCPU(dZeroEvtsRS, dZeroEvtsWS);
    timersub(&stopTime, &startTime, &totalTime);
    cpuTime += totalTime.tv_sec + totalTime.tv_usec / 1000000.0;
    fitRatioCPU(d0barEvtsRS, d0barEvtsWS);
    timersub(&stopTime, &startTime, &totalTime);
    cpuTime += totalTime.tv_sec + totalTime.tv_usec / 1000000.0;

    std::cout << "GPU time [seconds] : " << gpuTime << "\nCPU time [seconds] : " << cpuTime << std::endl;

    return 0;
}
