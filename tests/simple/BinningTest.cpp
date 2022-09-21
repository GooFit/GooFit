#include <goofit/Catch.h>

#include <goofit/BinnedDataSet.h>
#include <goofit/FitManager.h>
#include <goofit/PDFs/basic/ExpPdf.h>
#include <goofit/PDFs/basic/GaussianPdf.h>
#include <goofit/PDFs/combine/ProdPdf.h>
#include <goofit/UnbinnedDataSet.h>

#include <goofit/Variable.h>

#include <random>

using namespace GooFit;

TEST_CASE("Variations on a exp fit", "[simple][fit][exp]") {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Observable xvar{"xvar", 0, 10};

    std::unique_ptr<DataSet> data;
    size_t N;

    // Data set
    SECTION("Binned") {
        data.reset(new BinnedDataSet(xvar));
        N = 100000;
    }
    SECTION("Unbinned") {
        data.reset(new UnbinnedDataSet(xvar));
        N = 10000;
    }

    // Generate toy events.
    for(int i = 0; i < N; ++i) {
        double val = d(gen);
        if(val < 10) {
            xvar.setValue(val);
            data->addEvent();
        }
    }

    // Fit parameter
    Variable alpha{"alpha", -2, 0.1, -10, 10};

    // GooPdf object
    ExpPdf exppdf{"exppdf", xvar, alpha};
    exppdf.setData(data.get());

    FitManager fitter{&exppdf};
    fitter.fit();
    CHECK(fitter);

    double err = 1 / sqrt(N) * fabs(alpha.getValue());
    CHECK(alpha.getError() < 3 * err);
    CHECK(alpha.getValue() == Approx(-1.5).epsilon(3 / sqrt(N)));
}

TEST_CASE("Variations on a dual fit", "[simple][fit][exp]") {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> dx(1.5);
    std::exponential_distribution<> dy(.75);

    // Independent variable.
    Observable xvar{"xvar", 0, 10};
    Observable yvar{"yvar", 0, 10};

    // Data set
    std::unique_ptr<DataSet> data;
    size_t N;

    // Data set
    SECTION("Binned") {
        data.reset(new BinnedDataSet({xvar, yvar}));
        N = 200000;
    }
    SECTION("Unbinned") {
        data.reset(new UnbinnedDataSet({xvar, yvar}));
        N = 20000;
    }

    // Generate toy events.
    for(int i = 0; i < N; ++i) {
        double xval = dx(gen);
        double yval = dy(gen);
        if(xval < 10 && yval < 10) {
            xvar.setValue(xval);
            yvar.setValue(yval);
            data->addEvent();
        }
    }

    // Fit parameter
    Variable xalpha{"xalpha", -2, 0.1, -10, 10};
    // Fit parameter
    Variable yalpha{"yalpha", -2, 0.1, -10, 10};

    // GooPdf object
    ExpPdf xpdf{"xpdf", xvar, xalpha};
    ExpPdf ypdf{"ypdf", yvar, yalpha};
    ProdPdf totalpdf{"totalpdf", {&xpdf, &ypdf}};
    totalpdf.setData(data.get());

    FitManager fitter{&totalpdf};
    fitter.fit();

    CHECK(fitter);
    double xerr = 1 / sqrt(N) * fabs(xalpha.getValue());
    double yerr = 1 / sqrt(N) * fabs(yalpha.getValue());
    CHECK(xalpha.getError() < 3. * xerr);
    CHECK(yalpha.getError() < 3. * yerr);
    CHECK(xalpha.getValue() == Approx(-1.5).epsilon(3 / sqrt(N)));
    CHECK(yalpha.getValue() == Approx(-.75).epsilon(3 / sqrt(N)));
}

TEST_CASE("Dual fit with different fitter variable", "[simple][fit][gauss]") {
    // Random number generation
    std::mt19937 gen(137);
    std::normal_distribution<> dx(1.5, .3);
    std::normal_distribution<> dy(.75, .5);

    // Independent variable.
    Observable xvar{"xvar", -10, 10};
    Observable yvar{"yvar", -10, 10};
    //
    // Data set
    std::unique_ptr<DataSet> data;
    size_t N;
    bool binned;

    // Data set
    SECTION("Binned") {
        data.reset(new BinnedDataSet({xvar, yvar}, "Some name"));
        N      = 200000;
        binned = true;
    }
    SECTION("Unbinned") {
        data.reset(new UnbinnedDataSet({xvar, yvar}, "Some name"));
        N      = 20000;
        binned = false;
    }

    // Generate toy events.
    for(int i = 0; i < N; ++i) {
        xvar = dx(gen);
        yvar = dy(gen);
        if(xvar && yvar)
            data->addEvent();
    }

    // Fit parameter
    Variable xalpha{"xalpha", 0, 0.1, -10, 10};
    Variable xsigma{"xsigma", .1, 0.1, 0, 3};
    Variable unused{"unused", 1, .1, 0, 2};
    // Fit parameter
    Variable ysigma{"ysigma", .1, 0.1, 0, 3};
    Variable yalpha{"yalpha", 0, 0.1, -10, 10};

    // GooPdf object
    GaussianPdf ypdf{"ypdf", yvar, yalpha, ysigma};
    GaussianPdf xpdf{"xpdf", xvar, xalpha, xsigma};
    ProdPdf totalpdf{"totalpdf", {&xpdf, &ypdf}};
    totalpdf.setData(data.get());

    FitManager fitter{&totalpdf};
    fitter.fit();

    CHECK(fitter);
    double xaerr = 1 / sqrt(N) * fabs(xalpha.getValue());
    double yaerr = 1 / sqrt(N) * fabs(yalpha.getValue());
    double xserr = 1 / sqrt(N) * fabs(xsigma.getValue());
    double yserr = 1 / sqrt(N) * fabs(ysigma.getValue());
    CHECK(xalpha.getError() < 3. * xaerr);
    CHECK(yalpha.getError() < 3. * yaerr);
    CHECK(xsigma.getError() < 3. * xserr);
    CHECK(ysigma.getError() < 3. * yserr);

    CAPTURE(1 / sqrt(N));
    CAPTURE(binned);

    CHECK(xalpha.getValue() == Approx(1.5).epsilon(3 / sqrt(N)));
    CHECK(yalpha.getValue() == Approx(.75).epsilon(3 / sqrt(N)));
    CHECK(xsigma.getValue() == Approx(.3).epsilon(3 / sqrt(N) * (binned ? 10 : 1)));
    CHECK(ysigma.getValue() == Approx(.5).epsilon(3 / sqrt(N) * (binned ? 10 : 1)));
}
