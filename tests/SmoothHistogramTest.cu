#include "gtest/gtest.h"

#include "goofit/fitting/FitManagerMinuit1.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/PDFs/basic/SmoothHistogramPdf.h"

#include "goofit/Variable.h"

#include <sys/time.h>
#include <sys/times.h>
#include <iostream>

#include <random>

using namespace std;
using namespace GooFit;

TEST(SmoothHistogram, SimpleFit) {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Variable xvar{"xvar", 0, 10};

    // Data set
    BinnedDataSet data(&xvar);

    // Generate toy events.
    for(int i = 0; i < 10000; ++i) {
        double val = d(gen);
        if(val < 10) {
            xvar.setValue(val);
            data.addEvent();
        }
    }

    // Fit parameter
    Variable smoothing{"smooth", 1, 0, 3};

    // GooPdf object
    SmoothHistogramPdf pdf{"smoothhistogram", &data, &smoothing};
    pdf.setData(&data);

    GooFit::FitManagerMinuit1 fitter{&pdf};
    fitter.setVerbosity(2);
    fitter.fit();

    EXPECT_TRUE(fitter);
    EXPECT_LT(smoothing.getError(), .9);
    EXPECT_NEAR(0.9, smoothing.getValue(), smoothing.getError() * 3);
}
