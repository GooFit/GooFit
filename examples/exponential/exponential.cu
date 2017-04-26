#include "goofit/Application.h"
#include "goofit/Variable.h"
#include "goofit/FitManager.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/PDFs/ExpPdf.h"
#include <CLI/Timer.hpp>
#include <iostream>

using namespace std;

int main(int argc, char** argv) {
    GooFit::Application app("Exponential example", argc, argv);

    try {
        app.run();
    } catch (const GooFit::ParseError &e) {
        return app.exit(e);
    }

    // Independent variable.
    Variable xvar{"xvar", 0, log(1 + RAND_MAX/2)};

    // Data set
    UnbinnedDataSet data(&xvar);

    // Generate toy events.
    
    CLI::Timer gen_timer{"Generating took"};
    for(int i=0; i<100000; ++i) {
        xvar.value = xvar.upperlimit - log(1+rand()/2);

        if(xvar.value >= 0)
            data.addEvent();
    }
    std::cout << rang::fg::blue << gen_timer << rang::style::reset << std::endl;

    // Fit parameter
    Variable alpha{"alpha", -2, 0.1, -10, 10};
    // GooPdf object
    ExpPdf exppdf{"exppdf", &xvar, &alpha};
    exppdf.setData(&data);

    FitManager fitter{&exppdf};
    fitter.fit();
    
    if(alpha.value < -1.01 || alpha.value > -0.99)
        return 1;

    return fitter;
}
