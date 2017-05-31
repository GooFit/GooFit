#include "goofit/Application.h"
#include "goofit/Variable.h"
#include "goofit/fitting/Params.h"
#include "goofit/fitting/FCN.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/PDFs/ExpPdf.h"

#include <Minuit2/MnUserParameters.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/FunctionMinimum.h>
#include "Minuit2/MnPrint.h"

#include <iostream>

using namespace std;
using namespace GooFit;

int main(int argc, char **argv) {
    GooFit::Application app("Exponential example", argc, argv);

    try {
        app.run();
    } catch(const GooFit::ParseError &e) {
        return app.exit(e);
    }

    Minuit2::MnPrint::SetLevel(3);

    // Independent variable.
    Variable xvar{"xvar", 0, log(1 + RAND_MAX / 2)};

    // Data set
    UnbinnedDataSet data(&xvar);

    // Generate toy events.
    for(int i = 0; i < 100000; ++i) {
        try {
            xvar.setValue(xvar.getUpperLimit() - log(1 + rand() / 2));
            data.addEvent();
        } catch(const GooFit::OutOfRange &) {
        }
    }

    // Fit parameter
    Variable alpha{"alpha", -2, 0.1, -10, 10};
    // GooPdf object
    ExpPdf exppdf{"exppdf", &xvar, &alpha};
    exppdf.setData(&data);

    GooFit::Params upar{exppdf};
    GooFit::FCN fcn{upar};

    Minuit2::MnMigrad migrad{fcn, upar};

    Minuit2::FunctionMinimum min = migrad();
    std::cout << "min= " << min << std::endl;

    return 0;
}
