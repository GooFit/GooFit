#include "goofit/Variable.h"
#include "goofit/FitManager.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/PDFs/BWPdf.h"
#include "goofit/PDFs/GaussianPdf.h"
#include "goofit/PDFs/ConvolutionPdf.h"
#include <iostream>

#include "TRandom.h"

using namespace std;

double cpu_bw(double x, double x0, double gamma) {
    double ret = gamma;
    ret /= (2*sqrt(M_PI));
    ret /= ((x-x0)*(x-x0) + 0.25*gamma*gamma);
    return ret;
}

int main(int argc, char** argv) {
    // Independent variable.
    Variable* xvar = new Variable("xvar", -10, 10);

    Variable* gamma = new Variable("gamma", 2, 0.1, 0.1, 5);
    Variable* sigma = new Variable("sigma", 1.5, 0.1, 0.1, 5);
    Variable* x0    = new Variable("x0", 0.2, 0.05, -1, 1);
    Variable* zero  = new Variable("zero", 0);

    TRandom donram(42);
    // Data set
    UnbinnedDataSet data(xvar);

    // Generate toy events.
    for(int i = 0; i < 100000; ++i) {
        xvar->value = donram.Uniform(20) - 10;

        double bwvalue = cpu_bw(xvar->value, x0->value, gamma->value);
        double roll = donram.Uniform() * (2.0 / (sqrt(M_PI)*gamma->value));

        if(roll > bwvalue) {
            --i;
            continue;
        }

        xvar->value += donram.Gaus(0, sigma->value);

        if((xvar->value < xvar->lowerlimit) || (xvar->value > xvar->upperlimit)) {
            --i;
            continue;
        }

        data.addEvent();
    }

    BWPdf* breit = new BWPdf("breit", xvar, x0, gamma);
    GaussianPdf* gauss = new GaussianPdf("gauss", xvar, zero, sigma);
    ConvolutionPdf* convolution = new ConvolutionPdf("convolution", xvar, breit, gauss);
    convolution->setData(&data);

    FitManager fitter(convolution);
    fitter.fit();

    return 0;
}
