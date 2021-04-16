#include <CLI/Timer.hpp>
#include <goofit/Application.h>
#include <goofit/FitManager.h>
#include <goofit/PDFs/basic/Legendre.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>
#include <iostream>
#include <vector>
#include <TH1F.h>
#include <TStyle.h>
#include <stdlib.h>
#include <TCanvas.h>

using namespace std;

void fitAndPlot(GooFit::LegendrePdf *total, GooFit::UnbinnedDataSet *data, TH1F &dataHist, GooFit::Observable xvar, const char *fname) {
    TCanvas *foo = 0;
    foo = new TCanvas();
    total->setData(data);
    GooFit::FitManager fitter(total);
    fitter.fit();
    if(!fitter)
        std::exit(fitter);

    TH1F pdfHist("pdfHist", "", xvar.getNumBins(), xvar.getLowerLimit(), xvar.getUpperLimit());
    pdfHist.SetStats(false);

    GooFit::UnbinnedDataSet grid = total->makeGrid();
    total->setData(&grid);
    std::vector<std::vector<double>> pdfVals = total->getCompProbsAtDataPoints();

    double totalPdf = 0;

    for(int i = 0; i < grid.getNumEvents(); ++i) {
        grid.loadEvent(i);
        pdfHist.Fill(xvar.getValue(), pdfVals[0][i]);
        totalPdf += pdfVals[0][i];
    }

    for(int i = 0; i < xvar.getNumBins(); ++i) {
        double val = pdfHist.GetBinContent(i + 1);
        val /= totalPdf;
        val *= data->getNumEvents();
        pdfHist.SetBinContent(i + 1, val);
    }

    // foo->SetLogy(true);
    dataHist.SetMarkerStyle(8);
    dataHist.SetMarkerSize(0.5);
    dataHist.Draw("e p");
    pdfHist.SetLineColor(kBlue);
    pdfHist.SetLineWidth(3);
    pdfHist.Draw("lsame");
    foo->SaveAs(fname);
}
int main(int argc, char **argv) {
	GooFit::Application app("Legendre example", argc, argv);

	GOOFIT_PARSE(app);

	// Independent variable.
	GooFit::Observable xvar{"xvar", -1, 1};

	//Histogram declarations 
	TH1F legHist("LegendreHist","",xvar.getNumBins(),xvar.getLowerLimit(),xvar.getUpperLimit());
	legHist.SetStats(false);

	// Data set
	GooFit::UnbinnedDataSet data(xvar);

	// Generate toy events
	CLI::Timer gen_timer{"Generating took"};
	for(int i = 0; i < 100000; ++i) {
		try {
			xvar.setValue(xvar.getLowerLimit() + (((double)rand()/RAND_MAX)*(xvar.getUpperLimit() - xvar.getLowerLimit())));
			data.addEvent();
			legHist.Fill(xvar.getValue());
		} catch(const GooFit::OutOfRange &) {
		}
	}
	std::cout << GooFit::magenta << gen_timer << GooFit::reset << std::endl;

	// Fit parameter
	GooFit::Variable weight0{"weight0",0,0.1,-1000,1000};
	GooFit::Variable weight1{"weight1",1,0.1,-1000,1000};
	GooFit::Variable weight2{"weight2",0,0.1,-1000,1000};
	weight0.setFixed(1);
	weight1.setFixed(1);
	weight2.setFixed(1);
	std::vector<GooFit::Variable> weights;
	weights.push_back(weight0);
	weights.push_back(weight1);
	weights.push_back(weight2);    
	unsigned int max = 3;
	// GooPdf object
	GooFit::LegendrePdf *Legpdf = new GooFit::LegendrePdf{"Legpdf", xvar, weights, max};
	Legpdf->setData(&data);
	fitAndPlot(Legpdf, &data, legHist, xvar, "legendre_fit_cpp.png");
//	GooFit::FitManager fitter{&Legpdf};
//	fitter.fit();
//The above 2 lines are commented out because they were previously used for fitting the legendre plot without saving it to an image, and they become redundant with the function fitAndPlot

//	return fitter;
	return 0;
}
