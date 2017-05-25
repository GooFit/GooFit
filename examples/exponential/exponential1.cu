#include "goofit/Application.h"
#include "goofit/Variable.h"
#include "goofit/fitting/FitManagerMinuit1.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/PDFs/ExpPdf.h"
#include <iostream>

using namespace std;
using namespace GooFit;

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
    for(int i=0; i<100000; ++i) {
        try {
            xvar.setValue( xvar.getUpperLimit() - log(1+rand()/2));
            data.addEvent();
        } catch (const GooFit::OutOfRange &) {}
    }

    // Fit parameter
    Variable alpha{"alpha", -2, 0.1, -10, 10};
    // GooPdf object
    ExpPdf exppdf{"exppdf", &xvar, &alpha};
    exppdf.setData(&data);

    GooFit::FitManagerMinuit1 fitter{&exppdf};
    fitter.fit();

    return 0;
}
