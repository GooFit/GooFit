#include <goofit/Catch.h>

#include "testhelpers.h"
#include <goofit/PDFs/basic/StepPdf.h>
#include <goofit/UnbinnedDataSet.h>

#include <goofit/Variable.h>

#include <random>

using namespace GooFit;

TEST_CASE("Step", "[convert][fit]") {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Observable xvar{"xvar", 1, 10};

    // Data set
    UnbinnedDataSet data(xvar);

    // Generate toy events.
    for(int i = 0; i < 1000; ++i) {
        double val = d(gen);
        if(val < 10 && val > 1) {
            xvar.setValue(val);
            data.addEvent();
        }
    }

    // Fit parameter
    Variable x0{"x0", 1, 0, 2};

    // GooPdf object
    StepPdf pdf{"steppdf", xvar, x0};
    pdf.setData(&data);

    bool fitter = test_fitter(&pdf);

    CHECK(fitter);
}
