#include <goofit/Catch.h>

#include <goofit/PDFs/basic/ExpPdf.h>
#include <goofit/PDFs/basic/GaussianPdf.h>
#include <goofit/PDFs/combine/ProdPdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/fitting/FitManagerMinuit1.h>

#include <goofit/Variable.h>

#include <random>

using namespace GooFit;

TEST_CASE("Minuit1 Exp fit", "[minuit1][root][fit][exp]") {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Observable xvar{"xvar", 0, 10};

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
    Variable alpha{"alpha", -2, 0.1, -10, 10};

    // GooPdf object
    ExpPdf exppdf{"exppdf", xvar, alpha};
    exppdf.setData(&data);

    GooFit::FitManagerMinuit1 fitter{&exppdf};
    fitter.setVerbosity(2);
    fitter.fit();

    CHECK(fitter);
    CHECK(alpha.getError() < .1);
    CHECK(alpha.getValue() == Approx(-1.5).margin(.3));
}

TEST_CASE("Minuit1 Dual fit", "[minuit1][root][fit][exp]") {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> dx(1.5);
    std::exponential_distribution<> dy(.75);

    // Independent variables.
    Observable xvar{"xvar", 0, 10};
    Observable yvar{"yvar", 0, 10};

    // Data set
    UnbinnedDataSet data{{xvar, yvar}};

    // Generate toy events.
    for(int i = 0; i < 20000; ++i) {
        double xval = dx(gen);
        double yval = dy(gen);
        if(xval < 10 && yval < 10) {
            xvar.setValue(xval);
            yvar.setValue(yval);
            data.addEvent();
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
    totalpdf.setData(&data);

    GooFit::FitManagerMinuit1 fitter{&totalpdf};
    fitter.setVerbosity(0);
    fitter.fit();

    CHECK(fitter);
    CHECK(xalpha.getError() < .1);
    CHECK(yalpha.getError() < .1);
    CHECK(xalpha.getValue() == Approx(-1.5).margin(.3));
    CHECK(yalpha.getValue() == Approx(-.75).margin(.3));
}

TEST_CASE("Minuit1 different fitter variable", "[minuit1][root][fit][exp]") {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> dx(1.5);
    std::exponential_distribution<> dy(.75);

    // Independent variable.
    Observable xvar{"xvar", 0, 10};
    Observable yvar{"yvar", 0, 10};

    // Data set
    UnbinnedDataSet data{{xvar, yvar}, "Some name"};

    // Generate toy events.
    for(int i = 0; i < 20000; ++i) {
        xvar = dx(gen);
        yvar = dy(gen);
        if(xvar && yvar) {
            data.addEvent();
        }
    }

    // Fit parameter
    Variable xalpha{"xalpha", -2, 0.1, -10, 10};
    // Fit parameter
    Variable yalpha{"yalpha", -2, 0.1, -10, 10};

    // GooPdf object
    ExpPdf ypdf{"ypdf", yvar, yalpha};
    ExpPdf xpdf{"xpdf", xvar, xalpha};
    ProdPdf totalpdf{"totalpdf", {&xpdf, &ypdf}};
    totalpdf.setData(&data);

    GooFit::FitManagerMinuit1 fitter{&totalpdf};
    fitter.setVerbosity(2);
    fitter.fit();

    CHECK(fitter);
    CHECK(xalpha.getError() < .1);
    CHECK(yalpha.getError() < .1);
    CHECK(xalpha.getValue() == Approx(-1.5).margin(.3));
    CHECK(yalpha.getValue() == Approx(-.65).margin(.3));
}

TEST_CASE("Minuit1 Using fitter constants", "[minuit1][root][fit][gauss]") {
    // Random number generation
    std::mt19937 gen(137);
    std::normal_distribution<> dx(1.5, .3);
    std::normal_distribution<> dy(-.75, .2);

    // Independent variable.
    Observable xvar{"xvar", -5, 5};
    Observable yvar{"yvar", -5, 5};

    // Data set
    UnbinnedDataSet data{{xvar, yvar}, "Some name"};

    // Generate toy events.
    for(int i = 0; i < 20000; ++i) {
        double xval = dx(gen);
        double yval = dy(gen);
        try {
            xvar.setValue(xval);
            yvar.setValue(yval);
            data.addEvent();
        } catch(const GooFit::OutOfRange &) {
        }
    }

    // Fit parameter
    Variable xalpha{"xalpha", 2, 0.1, -10, 10};
    Variable xsigma{"xsigma", .2};
    Variable yalpha{"yalpha", -1, 0.1, -10, 10};
    Variable ysigma{"ysigma", .3};

    // GooPdf object
    GaussianPdf xpdf{"xpdf", xvar, xalpha, xsigma};
    GaussianPdf ypdf{"ypdf", yvar, yalpha, ysigma};
    ProdPdf totalpdf{"totalpdf", {&xpdf, &ypdf}};
    totalpdf.setData(&data);

    GooFit::FitManagerMinuit1 fitter{&totalpdf};
    fitter.setVerbosity(2);
    fitter.fit();

    CHECK(fitter);
    CHECK(xalpha.getError() < .1);
    CHECK(yalpha.getError() < .1);
    CHECK(xalpha.getValue() == Approx(1.5).margin(.3));
    CHECK(yalpha.getValue() == Approx(-.75).margin(.3));
    CHECK(xsigma.getValue() == .2_a);
    CHECK(ysigma.getValue() == .3_a);
}
