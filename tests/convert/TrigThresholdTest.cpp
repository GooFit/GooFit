#include <goofit/Catch.h>

#include "testhelpers.h"
#include <goofit/PDFs/basic/TrigThresholdPdf.h>
#include <goofit/UnbinnedDataSet.h>

#include <goofit/Variable.h>

#include <random>

using namespace GooFit;

TEST_CASE("TrigThresholdUpper", "[convert][fit]") {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Observable xvar{"xvar", -10, 10};

    // Data set
    UnbinnedDataSet data(xvar);

    // Generate toy events.
    for(int i = 0; i < 1000; ++i) {
        double val = d(gen);
        if(val < 10) {
            xvar.setValue(val);
            data.addEvent();
        }
    }

    // Fit parameter
    Variable thres{"m", 1, 0.1, -10, 10};
    Variable trig{"s", 1, 0, 3};
    Variable linCons{"w", 1, 0, 3};

    // GooPdf object
    TrigThresholdPdf pdf{"trigthresholdpdf", xvar, thres, trig, linCons, true};
    pdf.setData(&data);

    bool fitter = test_fitter(&pdf);

    CHECK(fitter);
}

TEST_CASE("TrigThresholdLower", "[convert][fit]") {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Observable xvar{"xvar", -10, 10};

    // Data set
    UnbinnedDataSet data(xvar);

    // Generate toy events.
    for(int i = 0; i < 1000; ++i) {
        double val = d(gen);
        if(val < 10) {
            xvar.setValue(val);
            data.addEvent();
        }
    }

    // Fit parameter
    Variable thres{"m", 1, 0.1, -10, 10};
    Variable trig{"s", 1, 0, 3};
    Variable linCons{"w", 1, 0, 3};

    // GooPdf object
    TrigThresholdPdf pdf{"trigthresholdpdf", xvar, thres, trig, linCons, false};
    pdf.setData(&data);

    bool fitter = test_fitter(&pdf);

    CHECK(fitter);
}

TEST_CASE("TrigThresholdVeryUpper", "[convert][fit]") {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Observable xvar{"xvar", -10, 10};
    Observable yvar{"yvar", -10, 10};

    // Data set
    UnbinnedDataSet data{{xvar, yvar}};

    // Generate toy events.
    for(int i = 0; i < 100000; ++i) {
        double val = d(gen);
        if(val < 10) {
            xvar.setValue(val);
            yvar.setValue(val);
            data.addEvent();
        }
    }

    // Fit parameter
    Variable thres{"m", 1, 0.1, 1, 10};
    Variable trig{"s", 1, 0, 3};
    Variable linCons{"w", 1, 0, 3};
    Variable mass{"mass", 1, 0, 3};

    // GooPdf object
    TrigThresholdPdf pdf{"trigthresholdpdf", xvar, yvar, thres, trig, linCons, mass, true};
    pdf.setData(&data);

    bool fitter = test_fitter(&pdf);

    CHECK(fitter);
}

TEST_CASE("TrigThresholdVeryLower", "[convert][fit]") {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Observable xvar{"xvar", -10, 10};
    Observable yvar{"yvar", -10, 10};

    // Data set
    UnbinnedDataSet data{{xvar, yvar}};

    // Generate toy events.
    for(int i = 0; i < 100000; ++i) {
        double val = d(gen);
        if(val < 10) {
            xvar.setValue(val);
            yvar.setValue(val);
            data.addEvent();
        }
    }

    // Fit parameter
    Variable thres{"m", 1, 0.1, -10, 10};
    Variable trig{"s", 1, 0, 3};
    Variable linCons{"w", 1, 0, 3};
    Variable mass{"mass", 1, 0, 3};

    // GooPdf object
    TrigThresholdPdf pdf{"trigthresholdpdf", xvar, yvar, thres, trig, linCons, mass, false};
    pdf.setData(&data);

    bool fitter = test_fitter(&pdf);

    CHECK(fitter);
}
