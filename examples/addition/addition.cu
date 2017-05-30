#include "goofit/Application.h"
#include "goofit/Variable.h"
#include "goofit/PDFs/GaussianPdf.h"
#include "goofit/PDFs/AddPdf.h"
#include "goofit/PDFs/PolynomialPdf.h"
#include "goofit/FitManager.h"
#include "goofit/UnbinnedDataSet.h"

#include "TRandom.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TCanvas.h"

#include <sys/time.h>
#include <sys/times.h>
#include <iostream>

using namespace std;
using namespace GooFit;

int main(int argc, char** argv) {
    GooFit::Application app("Addition example", argc, argv);

    try {
        app.run();
    } catch (const GooFit::ParseError &e) {
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

    vector<Variable*> vars;
    Variable* xvar = new Variable("xvar", -5, 5);
    vars.push_back(xvar);
    UnbinnedDataSet data(vars);

    TH1F xvarHist("xvarHist", "",
                  xvar->getNumBins(), xvar->getLowerLimit(), xvar->getUpperLimit());

    xvarHist.SetStats(false);

    TRandom donram(42);
    double totalData = 0;

    for(int i = 0; i < 100000; ++i) {
        xvar->setValue(donram.Gaus(0.2, 1.1));

        if(donram.Uniform() < 0.1)
            xvar->setValue(donram.Uniform(xvar->getLowerLimit(), xvar->getUpperLimit()));

        if(fabs(xvar->getValue()) > 5) {
            --i;
            continue;
        }

        data.addEvent();
        xvarHist.Fill(xvar->getValue());
        totalData++;
    }

    Variable* xmean = new Variable("xmean", 0, 1, -10, 10);
    Variable* xsigm = new Variable("xsigm", 1, 0.5, 1.5);
    GaussianPdf signal("signal", xvar, xmean, xsigm);

    vars.clear();
    Variable* constant = new Variable("constant", 1.0);
    vars.push_back(constant);
    PolynomialPdf backgr("backgr", xvar, vars);

    vector<PdfBase*> comps;
    comps.push_back(&signal);
    comps.push_back(&backgr);

    vars.clear();
    Variable* sigFrac = new Variable("sigFrac", 0.9, 0.75, 1.00);
    vars.push_back(sigFrac);

    AddPdf total("total", vars, comps);
    total.setData(&data);
    FitManager fitter(&total);
    fitter.fit();

    TH1F pdfHist("pdfHist", "",
                 xvar->getNumBins(), xvar->getLowerLimit(), xvar->getUpperLimit());
    TH1F sigHist("sigHist", "",
                 xvar->getNumBins(), xvar->getLowerLimit(), xvar->getUpperLimit());
    TH1F bkgHist("bkgHist", "",
                 xvar->getNumBins(), xvar->getLowerLimit(), xvar->getUpperLimit());

    pdfHist.SetStats(false);
    sigHist.SetStats(false);
    bkgHist.SetStats(false);

    UnbinnedDataSet grid(xvar);

    for(int i = 0; i < xvar->getNumBins(); ++i) {
        double step = (xvar->getUpperLimit() - xvar->getLowerLimit())/xvar->getNumBins();
        xvar->setValue(xvar->getLowerLimit() + (i + 0.5) * step);
        grid.addEvent();
    }

    total.setData(&grid);
    std::vector<std::vector<double>> pdfVals = total.getCompProbsAtDataPoints();

    TCanvas foo;

    double totalPdf = 0;

    for(int i = 0; i < grid.getNumEvents(); ++i) {
        grid.loadEvent(i);
        pdfHist.Fill(xvar->getValue(), pdfVals[0][i]);
        sigHist.Fill(xvar->getValue(), pdfVals[1][i]);
        bkgHist.Fill(xvar->getValue(), pdfVals[2][i]);
        totalPdf += pdfVals[0][i];
    }

    for(int i = 0; i < xvar->getNumBins(); ++i) {
        double val = pdfHist.GetBinContent(i+1);
        val /= totalPdf;
        val *= totalData;
        pdfHist.SetBinContent(i+1, val);
        val = sigHist.GetBinContent(i+1);
        val /= totalPdf;
        val *= sigFrac->getValue();
        val *= totalData;
        sigHist.SetBinContent(i+1, val);
        val = bkgHist.GetBinContent(i+1);
        val /= totalPdf;
        val *= (1.0 - sigFrac->getValue());
        val *= totalData;
        bkgHist.SetBinContent(i+1, val);
    }

    xvarHist.SetMarkerStyle(8);
    xvarHist.SetMarkerSize(0.5);
    xvarHist.Draw("p");
    pdfHist.SetLineColor(kBlue);
    pdfHist.SetLineWidth(3);
    pdfHist.Draw("lsame");
    sigHist.SetLineColor(kBlue);
    sigHist.SetLineStyle(kDashed);
    sigHist.SetLineWidth(3);
    sigHist.Draw("lsame");
    bkgHist.SetLineColor(kRed);
    bkgHist.SetLineWidth(3);
    bkgHist.Draw("lsame");
    foo.SaveAs("xhist.png");



    return 0;
}
