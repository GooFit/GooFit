#include <CLI/Timer.hpp>
#include <goofit/Application.h>
#include <goofit/FitManager.h>
#include <goofit/PDFs/basic/chebyshevPdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>
#include <iostream>
#include <vector>
#include <TH1F.h>
#include <TStyle.h>
#include <stdlib.h>
#include <TCanvas.h>

using namespace std;

void fitAndPlot(GooFit::chebyshevPdf *total,
                GooFit::UnbinnedDataSet *data,
                TH1F &dataHist,
                GooFit::Observable xvar,
                const char *fname) {
    TCanvas *foo = 0;
    foo          = new TCanvas();
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
    pdfHist.Draw("same");
    foo->SaveAs(fname);
}
int main(int argc, char **argv) {
	GooFit::Application app("chebyshevPdf example", argc, argv);

	GOOFIT_PARSE(app);
	
	// Independent variable.
	GooFit::Observable xvar{"xvar", -1, 1};
	xvar.setNumBins(25);

	// Histogram declarations
	double x[100000];
	double f, fmax = 0.0;
	TH1F chebHist("ChebyshevHist","",xvar.getNumBins(),xvar.getLowerLimit(),xvar.getUpperLimit());
	chebHist.SetStats(false);
	// Data set
	GooFit::UnbinnedDataSet data(xvar);

	// Generate toy events
	CLI::Timer gen_timer{"Generating took"};
	for(int i = 0; i < 100000; i++) {
		f = (double)rand()/RAND_MAX;
		x[i] = xvar.getLowerLimit()+ f*(xvar.getUpperLimit()-xvar.getLowerLimit());
		if(fmax < 1 + (2*x[i]*x[i])) {
			fmax = 1 + (2*x[i]*x[i]);
		}
	}
	for(int k = 0; k < 100000; k++) {
		f = (double)rand()/RAND_MAX;
		if(f < (1+(2*x[k]*x[k]))/fmax) {
			xvar.setValue(x[k]);
			data.addEvent();
			chebHist.Fill(xvar.getValue());
		}
	}
	std::cout << GooFit::magenta << gen_timer << GooFit::reset << std::endl;

	// Fit parameter
	GooFit::Variable weight0{"weight0",2,0.1,-1000,1000};
	GooFit::Variable weight1{"weight1",0,0.1,-1000,1000};
	GooFit::Variable weight2{"weight2",1.5,0.1,-1000,1000};
	GooFit::Variable weight3{"weight3",0,0.1,-1000,1000};
	weight0.setFixed(1);
	weight1.setFixed(1);
	weight2.setFixed(0);
	weight3.setFixed(1);
	std::vector<GooFit::Variable> weights;
	weights.push_back(weight0);
	weights.push_back(weight1);
	weights.push_back(weight2);    
	weights.push_back(weight3);    
	unsigned int max = 3;
	// Example data is generated as T_2 + 2*T_0. Equivalent to 2x^2 + 1
	// GooPdf object
	GooFit::chebyshevPdf *chebpdf = new GooFit::chebyshevPdf{"chebpdf", xvar, weights, max};
	chebpdf->setData(&data);
	fitAndPlot(chebpdf, &data, chebHist, xvar, "chebyshev_fit_cpp.png");
//	GooFit::FitManager fitter{&chebpdf};
//	fitter.fit();
	return 0;

}
