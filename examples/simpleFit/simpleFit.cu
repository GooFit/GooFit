#include "goofit/Application.h"
#include "goofit/FitManager.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/PDFs/LandauPdf.h"
#include "goofit/PDFs/NovosibirskPdf.h"
#include "goofit/PDFs/BifurGaussPdf.h"

#include "goofit/Variable.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TRandom.h"

#include <sys/time.h>
#include <sys/times.h>
#include <iostream>

using namespace std;
using namespace GooFit;

// CPU-side Novosibirsk evaluation for use in generating toy MC.
double novosib(double x, double peak, double width, double tail) {
    double qa = 0, qb = 0, qc = 0, qx = 0, qy = 0;

    if(fabs(tail) < 1.e-7)
        qc = 0.5 * pow(((x - peak) / width), 2);
    else {
        qa = tail * sqrt(log(4.));
        qb = sinh(qa) / qa;
        qx = (x - peak) / width * qb;
        qy = 1. + tail * qx;

        //---- Cutting curve from right side

        if(qy > 1.E-7)
            qc = 0.5 * (pow((log(qy) / tail), 2) + tail * tail);
        else
            qc = 15.0;
    }

    //---- Normalize the result

    return exp(-qc);
}

TCanvas *foo = 0;

void fitAndPlot(GooPdf *total, UnbinnedDataSet *data, TH1F &dataHist, Variable *xvar, const char *fname) {
    total->setData(data);
    FitManager fitter(total);
    fitter.fit();

    if(!fitter)
        std::exit(fitter);

    TH1F pdfHist("pdfHist", "", xvar->getNumBins(), xvar->getLowerLimit(), xvar->getUpperLimit());
    pdfHist.SetStats(false);

    UnbinnedDataSet grid(xvar);
    double step = (xvar->getUpperLimit() - xvar->getLowerLimit()) / xvar->getNumBins();

    for(int i = 0; i < xvar->getNumBins(); ++i) {
        xvar->setValue(xvar->getLowerLimit() + (i + 0.5) * step);
        grid.addEvent();
    }

    total->setData(&grid);
    std::vector<std::vector<double>> pdfVals = total->getCompProbsAtDataPoints();

    double totalPdf = 0;

    for(int i = 0; i < grid.getNumEvents(); ++i) {
        grid.loadEvent(i);
        pdfHist.Fill(xvar->getValue(), pdfVals[0][i]);
        totalPdf += pdfVals[0][i];
    }

    for(int i = 0; i < xvar->getNumBins(); ++i) {
        double val = pdfHist.GetBinContent(i + 1);
        val /= totalPdf;
        val *= data->getNumEvents();
        pdfHist.SetBinContent(i + 1, val);
    }

    // foo->SetLogy(true);
    dataHist.SetMarkerStyle(8);
    dataHist.SetMarkerSize(0.5);
    dataHist.Draw("p");
    pdfHist.SetLineColor(kBlue);
    pdfHist.SetLineWidth(3);
    pdfHist.Draw("lsame");
    foo->SaveAs(fname);
}

int main(int argc, char **argv) {
    GooFit::Application app("Simple fit example", argc, argv);

    size_t numevents = 100000;
    app.add_option("-n,--num", numevents, "Number of events", true);

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

    // Independent variable.
    Variable *xvar = new Variable("xvar", -100, 100);
    xvar->setNumBins(1000); // For such a large range, want more bins for better accuracy in normalisation.

    // Data sets for the three fits.
    UnbinnedDataSet landdata(xvar);
    UnbinnedDataSet bifgdata(xvar);
    UnbinnedDataSet novodata(xvar);

    // Histograms for showing the fit.
    TH1F landHist("landHist", "", xvar->getNumBins(), xvar->getLowerLimit(), xvar->getUpperLimit());
    TH1F bifgHist("bifgHist", "", xvar->getNumBins(), xvar->getLowerLimit(), xvar->getUpperLimit());
    TH1F novoHist("novoHist", "", xvar->getNumBins(), xvar->getLowerLimit(), xvar->getUpperLimit());
    landHist.SetStats(false);
    bifgHist.SetStats(false);
    novoHist.SetStats(false);

    TRandom donram(42);

    double maxNovo = 0;

    for(double x = xvar->getLowerLimit(); x < xvar->getUpperLimit(); x += 0.01) {
        double curr = novosib(x, 0.3, 0.5, 1.0);

        if(curr < maxNovo)
            continue;

        maxNovo = curr;
    }

    double leftSigma     = 13;
    double rightSigma    = 29;
    double leftIntegral  = 0.5 / (leftSigma * sqrt(2 * M_PI));
    double rightIntegral = 0.5 / (rightSigma * sqrt(2 * M_PI));
    double totalIntegral = leftIntegral + rightIntegral;
    double bifpoint      = -10;

    // Generating three sets of toy MC.
    while(landdata.getNumEvents() < numevents) {
        // Landau
        try {
            xvar->setValue(donram.Landau(20, 1));
            landdata.addEvent();
            landHist.Fill(xvar->getValue());
        } catch(const GooFit::OutOfRange &) {
        }
    }

    while(bifgdata.getNumEvents() < numevents) {
        // Bifurcated Gaussian
        double val;
        if(donram.Uniform() < (leftIntegral / totalIntegral)) {
            do {
                val = donram.Gaus(bifpoint, rightSigma);
            } while(val < bifpoint || val > xvar->getUpperLimit());
            xvar->setValue(val);

        } else {
            do {
                val = donram.Gaus(bifpoint, leftSigma);
            } while(val > bifpoint || val < xvar->getLowerLimit());
            xvar->setValue(val);
        }

        bifgdata.addEvent();
        bifgHist.Fill(xvar->getValue());
    }

    while(novodata.getNumEvents() < numevents) {
        // And Novosibirsk.
        while(true) {
            xvar->setValue(donram.Uniform(xvar->getLowerLimit(), xvar->getUpperLimit()));
            double y = donram.Uniform(0, maxNovo);

            if(y < novosib(xvar->getValue(), 0.3, 0.5, 1.0))
                break;
        }

        novodata.addEvent();
        novoHist.Fill(xvar->getValue());
    }

    foo = new TCanvas();

    Variable *mpv   = new Variable("mpv", 40, 0, 150);
    Variable *sigma = new Variable("sigma", 5, 0, 30);
    GooPdf *landau  = new LandauPdf("landau", xvar, mpv, sigma);
    fitAndPlot(landau, &landdata, landHist, xvar, "landau.png");

    Variable *nmean = new Variable("nmean", 0.4, -10.0, 10.0);
    Variable *nsigm = new Variable("nsigm", 0.6, 0.0, 1.0);
    Variable *ntail = new Variable("ntail", 1.1, 0.1, 0.0, 3.0);
    GooPdf *novo    = new NovosibirskPdf("novo", xvar, nmean, nsigm, ntail);
    fitAndPlot(novo, &novodata, novoHist, xvar, "novo.png");

    Variable *gmean = new Variable("gmean", 3.0, 1, -15, 15);
    Variable *lsigm = new Variable("lsigm", 10, 1, 10, 20);
    Variable *rsigm = new Variable("rsigm", 20, 1, 10, 40);
    GooPdf *bifur   = new BifurGaussPdf("bifur", xvar, gmean, lsigm, rsigm);
    fitAndPlot(bifur, &bifgdata, bifgHist, xvar, "bifur.png");

    return 0;
}
