#include <gtest/gtest.h>

#include "testhelpers.h"
#include <goofit/PDFs/basic/NovosibirskPdf.h>
#include <goofit/UnbinnedDataSet.h>

#include <goofit/Variable.h>

#include <random>

using namespace std;
using namespace GooFit;

TEST(Novosibirsk, SimpleFit) {
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

    EXPECT_TRUE(fitter);
    EXPECT_LT(m.getError(), .1);
    EXPECT_NEAR(0.1453, m.getValue(), m.getError() * 3);
    EXPECT_LT(s.getError(), .1);
    EXPECT_NEAR(0.2650, s.getValue(), s.getError() * 3);
    EXPECT_LT(t.getError(), .1);
    EXPECT_NEAR(0.8626, t.getValue(), t.getError() * 3);
}
