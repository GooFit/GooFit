#include <goofit/Catch.h>

#include "testhelpers.h"
#include <goofit/PDFs/basic/NovosibirskPdf.h>
#include <goofit/UnbinnedDataSet.h>

#include <goofit/Variable.h>

#include <random>

using namespace GooFit;

TEST_CASE("Novosibirsk", "[convert][fit]") {
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
    Variable m{"m", 1, 0.1, -10, 10};
    Variable s{"s", 1, 0, 3};
    Variable t{"t", 1, 0, 3};

    // GooPdf object
    NovosibirskPdf novopdf{"novopdf", xvar, m, s, t};
    novopdf.setData(&data);

    bool fitter = test_fitter(&novopdf);

    CHECK(fitter);
    CHECK(m.getError() < .1);
    CHECK(m.getValue() == Approx(0.1453).margin(m.getError() * 3));
    CHECK(s.getError() < .1);
    CHECK(s.getValue() == Approx(0.2650).margin(s.getError() * 3));
    CHECK(t.getError() < .1);
    CHECK(t.getValue() == Approx(0.8626).margin(t.getError() * 3));
}
