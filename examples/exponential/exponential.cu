#include "goofit/Application.h"
#include "goofit/Variable.h"
#include "goofit/FitManager.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/PDFs/ExpPdf.h"
#include <iostream>

using namespace std;

int main(int argc, char** argv) {
    GooFit::Application app("Exponential example", argc, argv);

    try {
        app.run();
    } catch (const GooFit::Error &e) {
        return app.exit(e);
    }

    // Independent variable.
    Variable* xvar = new Variable("xvar", 0, log(1 + RAND_MAX/2));

    // Data set
    UnbinnedDataSet data(xvar);

    // Generate toy events.
    for(int i = 0; i < 100000; ++i) {
        xvar->value = xvar->upperlimit - log(1+rand()/2);

        if(xvar->value < 0)
            continue;

        data.addEvent();
    }

    // Fit parameter
    Variable* alpha = new Variable("alpha", -2, 0.1, -10, 10);
    // GooPdf object
    ExpPdf* exppdf = new ExpPdf("exppdf", xvar, alpha);
    exppdf->setData(&data);

    FitManager fitter(exppdf);
    fitter.fit();

    return 0;
}
