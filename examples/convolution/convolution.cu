#include <goofit/Application.h>
#include <goofit/FitManager.h>
#include <goofit/PDFs/basic/BWPdf.h>
#include <goofit/PDFs/basic/GaussianPdf.h>
#include <goofit/PDFs/combine/ConvolutionPdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TRandom.h>

using namespace std;
using namespace GooFit;

double cpu_bw(double x, double x0, double gamma) {
    double ret = gamma;
    ret /= (2 * sqrt(M_PI));
    ret /= ((x - x0) * (x - x0) + 0.25 * gamma * gamma);
    return ret;
}

int main(int argc, char **argv) {
    GooFit::Application app("Convolution example", argc, argv);

    GOOFIT_PARSE(app);

    // Independent variable.
    Observable xvar{"xvar", -10, 10};

    Variable gamma{"gamma", 2, 0.1, 0.1, 5};
    Variable sigma{"sigma", 1.5, 0.1, 0.1, 5};
    Variable x0{"x0", 0.2, 0.05, -1, 1};
    Variable zero{"zero", 0};

    TRandom donram(42);
    // Data set
    UnbinnedDataSet data(xvar);

    // Generate toy events.
    for(int i = 0; i < 100000; ++i) {
        xvar.setValue(donram.Uniform(20) - 10);

        double bwvalue = cpu_bw(xvar.getValue(), x0.getValue(), gamma.getValue());
        double roll    = donram.Uniform() * (2.0 / (sqrt(M_PI) * gamma.getValue()));

        if(roll > bwvalue) {
            --i;
            continue;
        }

        xvar.setValue(xvar.getValue() + donram.Gaus(0, sigma.getValue()));

        if((xvar.getValue() < xvar.getLowerLimit()) || (xvar.getValue() > xvar.getUpperLimit())) {
            --i;
            continue;
        }

        data.addEvent();
    }

    BWPdf breit{"breit", xvar, x0, gamma};
    GaussianPdf gauss{"gauss", xvar, zero, sigma};
    ConvolutionPdf convolution{"convolution", xvar, &breit, &gauss};
    convolution.setData(&data);

    FitManager fitter(&convolution);
    fitter.fit();

    TFile f("output.root", "RECREATE");
    auto toroot = convolution.plotToROOT(xvar);
    toroot->Write();
    f.Write();

    return fitter;
}
