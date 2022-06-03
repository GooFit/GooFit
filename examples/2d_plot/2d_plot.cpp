#include <goofit/Application.h>
#include <goofit/FitManager.h>
#include <goofit/PDFs/basic/GaussianPdf.h>
#include <goofit/PDFs/combine/ProdPdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>
#include <goofit/utilities/Style.h>

#include <RVersion.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>

#include <iostream>
#include <random>
#include <sys/time.h>
#include <sys/times.h>

using namespace GooFit;

int main(int argc, char **argv) {
    GooFit::Application app("2D plot example", argc, argv);

    GOOFIT_PARSE(app);

    // In real code, use a random device here
    std::mt19937 gen(137);
    std::normal_distribution<> dx(0.2, 1.1);
    std::normal_distribution<> dy(0.5, 0.3);

    GooFit::setROOTStyle();

    Observable xvar{"xvar", -5, 5};
    Observable yvar{"yvar", -5, 5};
    UnbinnedDataSet data({xvar, yvar});

    TH2F dataHist("dataHist",
                  "",
                  xvar.getNumBins(),
                  xvar.getLowerLimit(),
                  xvar.getUpperLimit(),
                  yvar.getNumBins(),
                  yvar.getLowerLimit(),
                  yvar.getUpperLimit());
    TH1F xvarHist("xvarHist", "", xvar.getNumBins(), xvar.getLowerLimit(), xvar.getUpperLimit());
    TH1F yvarHist("yvarHist", "", yvar.getNumBins(), yvar.getLowerLimit(), yvar.getUpperLimit());

    dataHist.SetStats(false);
    xvarHist.SetStats(false);
    yvarHist.SetStats(false);

    size_t totalData = 0;

    while(totalData < 100000) {
        xvar.setValue(dx(gen));
        yvar.setValue(dy(gen));

        if(fabs(xvar.getValue()) < 5 && fabs(yvar.getValue()) < 5) {
            data.addEvent();
            dataHist.Fill(xvar.getValue(), yvar.getValue());
            xvarHist.Fill(xvar.getValue());
            yvarHist.Fill(yvar.getValue());
            totalData++;
        }
    }

    Variable xmean{"xmean", 0, 1, -10, 10};
    Variable xsigm{"xsigm", 1, 0.5, 1.5};
    GaussianPdf xgauss{"xgauss", xvar, xmean, xsigm};

    Variable ymean{"ymean", 0, 1, -10, 10};
    Variable ysigm{"ysigm", 0.4, 0.1, 0.6};
    GaussianPdf ygauss{"ygauss", yvar, ymean, ysigm};

    ProdPdf total("total", {&xgauss, &ygauss});
    total.setData(&data);

    FitManager fitter(&total);
    fitter.fit();

    TH2F pdfHist("pdfHist",
                 "",
                 xvar.getNumBins(),
                 xvar.getLowerLimit(),
                 xvar.getUpperLimit(),
                 yvar.getNumBins(),
                 yvar.getLowerLimit(),
                 yvar.getUpperLimit());
    TH1F xpdfHist("xpdfHist", "", xvar.getNumBins(), xvar.getLowerLimit(), xvar.getUpperLimit());
    TH1F ypdfHist("ypdfHist", "", yvar.getNumBins(), yvar.getLowerLimit(), yvar.getUpperLimit());

    pdfHist.SetStats(false);
    xpdfHist.SetStats(false);
    ypdfHist.SetStats(false);

    UnbinnedDataSet grid = total.makeGrid();
    total.setData(&grid);
    std::vector<std::vector<double>> pdfVals = total.getCompProbsAtDataPoints();

    TCanvas foo;
    dataHist.Draw("colz");
    foo.SaveAs("data.png");

    double totalPdf = 0;

    for(int i = 0; i < grid.getNumEvents(); ++i) {
        grid.loadEvent(i);
        pdfHist.Fill(xvar.getValue(), yvar.getValue(), pdfVals[0][i]);
        xpdfHist.Fill(xvar.getValue(), pdfVals[0][i]);
        ypdfHist.Fill(yvar.getValue(), pdfVals[0][i]);
        totalPdf += pdfVals[0][i];
    }

    for(int i = 0; i < xvar.getNumBins(); ++i) {
        double val = xpdfHist.GetBinContent(i + 1);
        val /= totalPdf;
        val *= totalData;
        xpdfHist.SetBinContent(i + 1, val);
    }

    for(int i = 0; i < yvar.getNumBins(); ++i) {
        double val = ypdfHist.GetBinContent(i + 1);
        val /= totalPdf;
        val *= totalData;
        ypdfHist.SetBinContent(i + 1, val);

        for(int j = 0; j < xvar.getNumBins(); ++j) {
            val = pdfHist.GetBinContent(j + 1, i + 1);
            val /= totalPdf;
            val *= totalData;
            pdfHist.SetBinContent(j + 1, i + 1, val);
        }
    }

    pdfHist.Draw("colz");
    foo.SaveAs("pdf.png");

    for(int i = 0; i < yvar.getNumBins(); ++i) {
        for(int j = 0; j < xvar.getNumBins(); ++j) {
            double pval = pdfHist.GetBinContent(j + 1, i + 1);
            double dval = dataHist.GetBinContent(j + 1, i + 1);
            pval -= dval;
            pval /= std::max(1.0, sqrt(dval));
            pdfHist.SetBinContent(j + 1, i + 1, pval);
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

    return fitter;
}
