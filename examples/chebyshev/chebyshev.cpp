
#include <CLI/Timer.hpp>
#include <goofit/Application.h>
#include <goofit/FitManager.h>
#include <goofit/PDFs/basic/chebyshevPdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>
#include <iostream>

int main(int argc, char **argv) {
	GooFit::Application app("chebyshevPdf example", argc, argv);

	GOOFIT_PARSE(app);

	// Independent variable.
	GooFit::Observable xvar{"xvar", -1, 1};

	// Data set
	GooFit::UnbinnedDataSet data(xvar);

	// Generate toy events
	CLI::Timer gen_timer{"Generating took"};
	for(int i = 0; i < 100000; ++i) {
		try {
			xvar.setValue(xvar.getLowerLimit() - (rand()*(xvar.getUpperLimit() - xvar.getLowerLimit())));
			data.addEvent();
		} catch(const GooFit::OutOfRange &) {
		}
	}

	std::cout << GooFit::magenta << gen_timer << GooFit::reset << std::endl;

	// Fit parameter
	GooFit::Variable weight0{"weight0",-2,0.1,-1000,1000};
	GooFit::Variable weight1{"weight1",-2,0.1,-1000,1000};
	GooFit::Variable weight2{"weight2",-2,0.1,-1000,1000};
	std::vector<GooFit::Variable> weights;
	weights.push_back(weight0);
	weights.push_back(weight1);
	weights.push_back(weight2);    
	unsigned int max = 3;
	// GooPdf object
	GooFit::chebyshevPdf chebpdf{"chebpdf", xvar, weights, max};
	chebpdf.setData(&data);

	GooFit::FitManager fitter{&chebpdf};
	fitter.fit();

	return fitter;
}
