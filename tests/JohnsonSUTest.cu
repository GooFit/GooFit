#include <gtest/gtest.h>

#include <goofit/PDFs/basic/JohnsonSUPdf.h>
#include <goofit/UnbinnedDataSet.h>
#include "testhelpers.h"

#include <goofit/Variable.h>





#include <random>

using namespace std;
using namespace GooFit;

TEST(JohnsonSU, SimpleFit) {
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
    Variable g{"g", 1, 0, 3};
    Variable e{"e", 1, 0, 3};

    // GooPdf object
    JohnsonSUPdf johnsonpdf{"johnsonpdf", xvar, m, s, g, e};
    johnsonpdf.setData(&data);

    bool fitter = test_fitter(&johnsonpdf);
    
    

    EXPECT_TRUE(fitter);
    EXPECT_LT(m.getError(), .1);
    EXPECT_NEAR(0.5, m.getValue(), m.getError() * 3);
}
