#include <CLI/Timer.hpp>
#include <goofit/Application.h>
#include <goofit/FitManager.h>
#include <goofit/PDFs/basic/poissonPdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <TCanvas.h>
#include <TH1F.h>
#include <TStyle.h>

using namespace std;


void fitAndPlot(GooFit::PoissonPdf *total, GooFit::UnbinnedDataSet *data, TH1F &dataHist, GooFit::Observable xvar, const char *fname) {
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
	GooFit::Application app("Poisson example", argc, argv);

	GOOFIT_PARSE(app);
	TCanvas *foo = 0;

	// Independent variable.
	GooFit::Observable xvar{"xvar", 0, 20}; //change xvar range? 
	xvar.setNumBins(20);
	TH1F poissonHist("PoissonHist","",xvar.getNumBins(),xvar.getLowerLimit(),xvar.getUpperLimit());
	poissonHist.SetStats(false);
	int ne = 0;
	double x[1000];
	const double e = 2.718281828459045235360287471352;
	vector<double> x_passed;
	const double xmin = 0.0;
	const double xmax = 20.0;
	double f = 0.0;
	double fmax = 0.0;
	double r = 0.0;
	const int Lambda = 10;
	//Data Set
	GooFit::UnbinnedDataSet data(xvar);
	
	// Generate toy events
	CLI::Timer gen_timer{"Generating took"};
	
	//Data generation
	for(int j = 0; j < 1000; j++) {
		f = (double)rand()/RAND_MAX;
		x[j] = xmin + f*(xmax-xmin);
		if (fmax < (pow(Lambda,x[j])*pow(e,-1*Lambda)/tgamma(x[j]+1))) { fmax=pow(Lambda,x[j])*pow(e,-1*Lambda)/tgamma(x[j]+1);}
	}
	std::cout << "fmax = " << fmax << endl;
	for (int k = 0; k < 1000; k++) { 
		r = (double)rand()/RAND_MAX;
		std::cout << "r = " << r << endl;
		std::cout << "vs - " << pow(Lambda,x[k])*pow(e,-1*Lambda)/(tgamma(x[k]+1)*fmax) << endl;
		if (r < (pow(Lambda,x[k])*pow(e,-1*Lambda)/(tgamma(x[k]+1)*fmax))) { 
			xvar.setValue(x[k]);
			data.addEvent();
			poissonHist.Fill(xvar.getValue());
			ne+=1;
		}
	}
	std::cout << "Number of events added: " << ne << endl;
	std::cout << GooFit::magenta << gen_timer << GooFit::reset << std::endl;
	// Fit parameter
	GooFit::Variable lambda{"lambda", 9, 0.001, 5, 15};
	// GooPdf object
	GooFit::PoissonPdf *poissonpdf = new GooFit::PoissonPdf{"poissonpdf", xvar, lambda};
	poissonpdf->setData(&data);
	fitAndPlot(poissonpdf, &data, poissonHist, xvar, "poisson_fit_cpp.png");
//	GooFit::FitManager fitter(poissonpdf);
//	fitter.fit();

	return 0; //Before commenting out the above 2 lines, it was return fitter;
}
