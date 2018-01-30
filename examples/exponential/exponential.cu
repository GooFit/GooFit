#include <CLI/Timer.hpp>
#include <goofit/Application.h>
#include <goofit/FitManager.h>
#include <goofit/PDFs/basic/ExpPdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>
#include <iostream>

int main(int argc, char **argv) {
    GooFit::Application app("Exponential example", argc, argv);

    GOOFIT_PARSE(app);

    // Independent variable.
    GooFit::Observable xvar{"xvar", 0, log(1 + RAND_MAX / 2)};

    // Data set
    GooFit::UnbinnedDataSet data(xvar);

    // Generate toy events
    CLI::Timer gen_timer{"Generating took"};
    for(int i = 0; i < 100000; ++i) {
        try {
            xvar.setValue(xvar.getUpperLimit() - log(1 + rand() / 2));
            data.addEvent();
        } catch(const GooFit::OutOfRange &) {
        }
    }

    std::cout << GooFit::magenta << gen_timer << GooFit::reset << std::endl;

    // Fit parameter
    GooFit::Variable alpha{"alpha", -2, 0.1, -10, 10};
    // GooPdf object
    GooFit::ExpPdf exppdf{"exppdf", xvar, alpha};
    exppdf.setData(&data);

    GooFit::FitManager fitter{&exppdf};
    fitter.fit();

    if(alpha.getValue() < -1.01 || alpha.getValue() > -0.99)
        return 1;

    return fitter;
}
