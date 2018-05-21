#include <goofit/Catch.h>

#include <random>

#include <goofit/FitManager.h>
#include <goofit/PDFs/basic/GaussianPdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>

#ifdef ROOT_FOUND
#include <goofit/fitting/FitManagerMinuit1.h>
#endif

using namespace GooFit;

TEST_CASE("No blind Gaussian fit", "[simple][gauss][fit]") {
    // Independent variable.
    Observable xvar{"xvar", 0, 10};

    // Data set
    UnbinnedDataSet data{xvar};

    // Random number generation
    std::mt19937 gen(137);
    std::normal_distribution<> d(1.5, .3);

    // Generate toy events.
    for(int i = 0; i < 100000; ++i) {
        double val = d(gen);
        if(std::fabs(val) < 10) {
            xvar.setValue(val);
            data.addEvent();
        }
    }

    // Fit parameter
    Variable alpha{"alpha", -2, -10, 10};
    Variable sigma{"sigma", 1, 0, 3};

    // GooPdf object
    GaussianPdf gausspdf{"gausspdf", xvar, alpha, sigma};
    gausspdf.fitTo(&data);

    CHECK(alpha.getError() < .01);
    CHECK(alpha.getValue() == Approx(1.5).margin(.03));
}

TEST_CASE("Blind Gaussian fit", "[blind][simple][gauss][fit]") {
    // Independent variable.
    Observable xvar{"xvar", 0, 10};

    // Data set
    UnbinnedDataSet data{xvar};

    // Random number generation
    std::mt19937 gen(137);
    std::normal_distribution<> d(1.5, .3);

    // Generate toy events.
    for(int i = 0; i < 100000; ++i) {
        double val = d(gen);
        if(std::fabs(val) < 10) {
            xvar.setValue(val);
            data.addEvent();
        }
    }

    // Fit parameter
    Variable alpha{"alpha", -2, -10, 10};
    Variable sigma{"sigma", 1, 0, 3};

    // Blinding (normally would be randomly generated)
    alpha.setBlind(1);

    // GooPdf object
    GaussianPdf gausspdf{"gausspdf", xvar, alpha, sigma};

    SECTION("Minuit 2") { gausspdf.fitTo(&data); }

#ifdef ROOT_FOUND
    SECTION("Minuit 1") {
        gausspdf.setData(&data);
        GooFit::FitManagerMinuit1 fitman{&gausspdf};
        fitman.fit();
        CHECK(fitman);
    }
#endif

    CHECK(alpha.getError() < .01);
    CHECK(alpha.getValue() == Approx(2.5).margin(.03));
}
