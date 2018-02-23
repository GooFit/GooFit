#include <gtest/gtest.h>

#include "testhelpers.h"
#include <goofit/PDFs/basic/SmoothHistogramPdf.h>
#include <goofit/UnbinnedDataSet.h>

#include <goofit/Variable.h>

#include <random>

using namespace std;
using namespace GooFit;

TEST(SmoothHistogram, SimpleFit) {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Observable xvar{"xvar", 0, 10};

    // Data set
    BinnedDataSet data(xvar);

    // Generate toy events.
    for(int i = 0; i < 10000; ++i) {
        double val = d(gen);
        if(val < 10 && val > 0) {
            xvar.setValue(val);
            data.addEvent();
        }
    }

    // Fit parameter
    Variable smoothing{"smooth", 1, 0, 1};

    // GooPdf object
    SmoothHistogramPdf pdf{"smoothhistogram", &data, smoothing};
    pdf.setData(&data);

    bool fitter = test_fitter(&pdf);

    EXPECT_TRUE(fitter);
    // EXPECT_LT(smoothing.getError(), .9);
    // EXPECT_NEAR(0.9, smoothing.getValue(), smoothing.getError() * 3);
}
