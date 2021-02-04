#include <CLI/Timer.hpp>
#include <goofit/Application.h>
#include <goofit/FitManager.h>
#include <goofit/PDFs/basic/poissonPdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <vector>

using namespace std;

int main(int argc, char **argv) {
	GooFit::Application app("Poisson example", argc, argv);

	GOOFIT_PARSE(app);

	// Independent variable.
	GooFit::Observable xvar{"xvar", 0, 20}; //change xvar range? 
	double x[1000];
	const double e = 2.718281828459045235360287471352;
	std::vector<double> x_passed;
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
		if (fmax < (pow(Lambda,x[j])*pow(e,-1*Lambda)/tgamma(x[j]+1))) { fmax=x[j];}
	}
	for (int k = 0; k < 1000; k++) { 
		r = rand();
		if (r < (pow(Lambda,x[k])*pow(e,-1*Lambda)/(tgamma(x[k]+1)*fmax))) {
			xvar.setValue(x[k]);
			data.addEvent();
		}
	}

	std::cout << GooFit::magenta << gen_timer << GooFit::reset << std::endl;

	// Fit parameter
	GooFit::Variable lambda{"lambda", 3, 0.001, -11, 11};
	// GooPdf object
	GooFit::PoissonPdf poissonpdf{"poissonpdf", xvar, lambda};
	poissonpdf.setData(&data);

	GooFit::FitManager fitter{&poissonpdf};
	fitter.fit();

	return fitter;
}
