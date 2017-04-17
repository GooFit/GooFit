#include "goofit/Application.h"
#include "goofit/Variable.h"
#include "goofit/PDFs/GaussianPdf.h"
#include "goofit/PDFs/ProdPdf.h"
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

int main(int argc, char** argv) {

    GooFit::Application app("2D plot example", argc, argv);

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
    Variable* yvar = new Variable("yvar", -5, 5);
    vars.push_back(yvar);
    UnbinnedDataSet data(vars);

    TH2F dataHist("dataHist", "",
                  xvar->numbins, xvar->lowerlimit, xvar->upperlimit,
                  yvar->numbins, yvar->lowerlimit, yvar->upperlimit);
    TH1F xvarHist("xvarHist", "",
                  xvar->numbins, xvar->lowerlimit, xvar->upperlimit);
    TH1F yvarHist("yvarHist", "",
                  yvar->numbins, yvar->lowerlimit, yvar->upperlimit);

    dataHist.SetStats(false);
    xvarHist.SetStats(false);
    yvarHist.SetStats(false);

    TRandom donram(42);
    double totalData = 0;

    for(int i = 0; i < 100000; ++i) {
        xvar->value = donram.Gaus(0.2, 1.1);
        yvar->value = donram.Gaus(0.5, 0.3);

        if(fabs(xvar->value) > 5) {
            --i;
            continue;
        }

        if(fabs(yvar->value) > 5) {
            --i;
            continue;
        }

        data.addEvent();
        dataHist.Fill(xvar->value, yvar->value);
        xvarHist.Fill(xvar->value);
        yvarHist.Fill(yvar->value);
        totalData++;
    }

    Variable* xmean = new Variable("xmean", 0, 1, -10, 10);
    Variable* xsigm = new Variable("xsigm", 1, 0.5, 1.5);
    GaussianPdf xgauss("xgauss", xvar, xmean, xsigm);

    Variable* ymean = new Variable("ymean", 0, 1, -10, 10);
    Variable* ysigm = new Variable("ysigm", 0.4, 0.1, 0.6);
    GaussianPdf ygauss("ygauss", yvar, ymean, ysigm);

    vector<PdfBase*> comps;
    comps.push_back(&xgauss);
    comps.push_back(&ygauss);

    ProdPdf total("total", comps);
    total.setData(&data);
    FitManager fitter(&total);
    fitter.fit();

    TH2F pdfHist("pdfHist", "",
                 xvar->numbins, xvar->lowerlimit, xvar->upperlimit,
                 yvar->numbins, yvar->lowerlimit, yvar->upperlimit);
    TH1F xpdfHist("xpdfHist", "",
                  xvar->numbins, xvar->lowerlimit, xvar->upperlimit);
    TH1F ypdfHist("ypdfHist", "",
                  yvar->numbins, yvar->lowerlimit, yvar->upperlimit);

    pdfHist.SetStats(false);
    xpdfHist.SetStats(false);
    ypdfHist.SetStats(false);

    UnbinnedDataSet grid(vars);

    for(int i = 0; i < xvar->numbins; ++i) {
        double step = (xvar->upperlimit - xvar->lowerlimit)/xvar->numbins;
        xvar->value = xvar->lowerlimit + (i + 0.5) * step;

        for(int j = 0; j < yvar->numbins; ++j) {
            step = (yvar->upperlimit - yvar->lowerlimit)/yvar->numbins;
            yvar->value = yvar->lowerlimit + (j + 0.5) * step;
            grid.addEvent();
        }
    }

    total.setData(&grid);
    vector<vector<double>> pdfVals;
    total.getCompProbsAtDataPoints(pdfVals);

    TCanvas foo;
    dataHist.Draw("colz");
    foo.SaveAs("data.png");

    double totalPdf = 0;

    for(int i = 0; i < grid.getNumEvents(); ++i) {
        grid.loadEvent(i);
        pdfHist.Fill(xvar->value, yvar->value, pdfVals[0][i]);
        xpdfHist.Fill(xvar->value, pdfVals[0][i]);
        ypdfHist.Fill(yvar->value, pdfVals[0][i]);
        totalPdf += pdfVals[0][i];
    }

    for(int i = 0; i < xvar->numbins; ++i) {
        double val = xpdfHist.GetBinContent(i+1);
        val /= totalPdf;
        val *= totalData;
        xpdfHist.SetBinContent(i+1, val);
    }

    for(int i = 0; i < yvar->numbins; ++i) {
        double val = ypdfHist.GetBinContent(i+1);
        val /= totalPdf;
        val *= totalData;
        ypdfHist.SetBinContent(i+1, val);

        for(int j = 0; j < xvar->numbins; ++j) {
            val = pdfHist.GetBinContent(j+1, i+1);
            val /= totalPdf;
            val *= totalData;
            pdfHist.SetBinContent(j+1, i+1, val);
        }
    }

    pdfHist.Draw("colz");
    foo.SaveAs("pdf.png");

    for(int i = 0; i < yvar->numbins; ++i) {
        for(int j = 0; j < xvar->numbins; ++j) {
            double pval = pdfHist.GetBinContent(j+1, i+1);
            double dval = dataHist.GetBinContent(j+1, i+1);
            pval -= dval;
            pval /= std::max(1.0, sqrt(dval));
            pdfHist.SetBinContent(j+1, i+1, pval);
        }
    }

    pdfHist.GetZaxis()->SetRangeUser(-5, 5);
    pdfHist.Draw("colz");
    foo.SaveAs("pull.png");

    xvarHist.SetMarkerStyle(8);
    xvarHist.SetMarkerSize(0.5);
    xvarHist.Draw("p");
    xpdfHist.SetLineColor(kBlue);
    xpdfHist.SetLineWidth(3);
    xpdfHist.Draw("lsame");
    foo.SaveAs("xhist.png");

    yvarHist.SetMarkerStyle(8);
    yvarHist.SetMarkerSize(0.5);
    yvarHist.Draw("p");
    ypdfHist.SetLineColor(kBlue);
    ypdfHist.SetLineWidth(3);
    ypdfHist.Draw("lsame");
    foo.SaveAs("yhist.png");


    return 0;
}
