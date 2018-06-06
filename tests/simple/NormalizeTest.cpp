#include <goofit/Catch.h>

#include <goofit/FitManager.h>
#include <goofit/PDFs/basic/ExpPdf.h>
#include <goofit/PDFs/basic/GaussianPdf.h>
#include <goofit/PDFs/combine/ProdPdf.h>
#include <goofit/UnbinnedDataSet.h>

#include <goofit/Variable.h>

#include <random>

using namespace GooFit;

TEST_CASE("Known Normalize", "[basic]") {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Observable xvar{"xvar", 0, 10};

    // Data set
    UnbinnedDataSet data(xvar);

    // Generate toy events.
    for(int i = 0; i < 1000; ++i) {
        xvar = d(gen);
        if(xvar) {
            data.addEvent();
        }
    }

    // Fit parameter
    Variable alpha{"alpha", -2, 0.1, -10, 10};

    // GooPdf object
    ExpPdf exppdf{"exppdf", xvar, alpha};
    exppdf.setData(&data);
    exppdf.copyParams();

    CHECK(exppdf.integrate(0, 10) == Approx(0.5));
    CHECK(exppdf.normalize() == Approx(0.5));
    CHECK(exppdf.normalize() == Approx(0.5)); // Just verifying that it does not crash

    FitManager fitter{&exppdf};
    fitter.setVerbosity(0);
    fitter.fit();

    CHECK(alpha.getValue() == Approx(-1.5).epsilon(.05));
    CHECK(exppdf.normalize() == Approx(0.665099));
    CHECK(exppdf.normalize() == Approx(0.665099));
    // Just verifying that it does not crash
    CHECK(exppdf.normalize() == Approx(exppdf.integrate(0, 10)));
}

TEST_CASE("1D Grid", "[basic]") {
    const char *info =
        R"raw(Recreate in Python with:

    import numpy as np
    from scipy.stats import norm
    bin_centers = np.linspace(-5,10,101) + .15/2
    vv = norm.pdf(bin_centers,1.5,2.5)
    for i in range(0,80,20):
    print(i, bin_centers[i], vv[i]))raw";

    INFO(info);

    // Independent variable.
    Observable xvar{"xvar", -5, 10};

    // Fit parameter
    Variable mean{"mean", 1.5};
    Variable sigma{"sigma", 2.5};

    // GooPdf object
    GaussianPdf gausspdf{"gausspdf", xvar, mean, sigma};

    // Make a grid and use it
    auto dataset = gausspdf.makeGrid();
    gausspdf.setData(&dataset);

    CHECK(dataset.getValue(xvar, 0) == Approx(-4.925));
    CHECK(dataset.getValue(xvar, 20) == Approx(-1.925));
    CHECK(dataset.getValue(xvar, 40) == 1.075_a);
    CHECK(dataset.getValue(xvar, 60) == 4.075_a);

    // Compute probabilities at points
    auto vv = gausspdf.getCompProbsAtDataPoints();

    // Check
    REQUIRE(vv.size() == 1);
    REQUIRE(vv[0].size() == 100);

    CHECK(vv[0][0] == 0.00587129964482_a);
    CHECK(vv[0][20] == 0.0624318782242_a);
    CHECK(vv[0][40] == 0.157287605852_a);
    CHECK(vv[0][60] == 0.0938855055588_a);
}

/*
 */

TEST_CASE("2D Grid", "[basic]") {
    const char *info =
        R"raw(To recreate in Python:

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import multivariate_normal

    bin_centers_x = np.linspace(-5,10,101) + .15/2
    bin_centers_y = np.linspace(-5,10,101) + .15/2
    X,Y = np.meshgrid(bin_centers_x, bin_centers_y)

    vv = multivariate_normal.pdf(np.stack([X,Y],-1),[1.5, -1],[2.5, .5])

    for i in range(0,101,20):
    for j in range(0,101,20):
    print(i, j, X[i,j], Y[i,j], vv[i,j])

    plt.contourf(X,Y,vv,40))raw";

    INFO(info);

    // Independent variable.
    Observable xvar{"xvar", -5, 10};
    Observable yvar{"yvar", -5, 10};

    // Fit parameter
    Variable mean1{"mean1", 1.5};
    Variable sigma1{"sigma1", 2.5};
    Variable mean2{"mean2", -1};
    Variable sigma2{"sigma2", 0.5};

    // GooPdf object
    GaussianPdf gausspdf1{"gausspdf1", xvar, mean1, sigma1};
    GaussianPdf gausspdf2{"gausspdf2", yvar, mean2, sigma2};

    ProdPdf product("product", {&gausspdf1, &gausspdf2});

    // Make a grid and use it
    auto dataset = product.makeGrid();
    product.setData(&dataset);

    CHECK(dataset.getValue(xvar, 0) == Approx(-4.925));
    CHECK(dataset.getValue(xvar, 20) == Approx(-1.925));
    CHECK(dataset.getValue(xvar, 40) == 1.075_a);
    CHECK(dataset.getValue(xvar, 60) == 4.075_a);

    CHECK(dataset.getValue(yvar, 100 * 0) == Approx(-4.925));
    CHECK(dataset.getValue(yvar, 100 * 20) == Approx(-1.925));
    CHECK(dataset.getValue(yvar, 100 * 40) == 1.075_a);
    CHECK(dataset.getValue(yvar, 100 * 60) == 4.075_a);

    // Compute probabilities at points
    auto vv = product.getCompProbsAtDataPoints();

    // Check
    REQUIRE(vv.size() == 3);
    REQUIRE(vv[0].size() == 10000);
    REQUIRE(vv[1].size() == 10000);
    REQUIRE(vv[2].size() == 10000);

    CHECK(vv[1][0] == 0.00587129964482_a);
    CHECK(vv[1][20] == 0.0624318782242_a);
    CHECK(vv[1][40] == 0.157287605852_a);
    CHECK(vv[1][60] == 0.0938855055588_a);

    // EXPECT_FLOAT_EQ(vv[0][0], 7.53698529533e-12);
    // EXPECT_FLOAT_EQ(vv[0][20*100+20], 0.00579245632235);
    // EXPECT_FLOAT_EQ(vv[0][20*100+60], 0.0160636095436);
}

TEST_CASE("1D Eval", "[basic]") {
    // Independent variable.
    Observable xvar{"xvar", -5, 10};

    // Fit parameter
    Variable mean{"mean", 1.5};
    Variable sigma{"sigma", 2.5};

    // GooPdf object
    GaussianPdf gauss{"gausspdf", xvar, mean, sigma};

    auto v = gauss.evaluateAtPoints(xvar);

    CHECK(v[0] == 0.00587129964482_a);
    CHECK(v[20] == 0.0624318782242_a);
    CHECK(v[40] == 0.157287605852_a);
    CHECK(v[60] == 0.0938855055588_a);

    CHECK(gauss.normalize() == 6.26657_a);

    /*
    import numpy as np

    mean = 1.5
    sigma = 2.5
    x = np.array((1.2, 7.2))
    f = np.exp(-0.5*(x - mean)**2 / sigma**2)
    fn = f / np.sqrt(2 * np.pi * sigma**2)

    print("g({}) = {} -> {} normalized".format(x[0], f[0], fn[0]))
    print("g({}) = {} -> {} normalized".format(x[1], f[1], fn[1]))

    # g(1.2) = 0.9928258579038134 -> 0.158743208471746242 normalized
    # g(7.2) = 0.07433302085078965 -> 0.011861833938936505 normalized
    */

    xvar = 1.2;
    CHECK(gauss.getValue(EvalFunc::Prob) == 0.15843208471746245_a);
    CHECK(gauss.getValue(EvalFunc::Eval) == 0.99282585857903813_a);

    xvar = 7.2;
    CHECK(gauss.getValue(EvalFunc::Prob) == 0.0118618339389365_a);
    CHECK(gauss.getValue(EvalFunc::Eval) == 0.0743330208507896_a);
}

TEST_CASE("2D Eval", "[basic]") {
    // Independent variable.
    Observable xvar{"xvar", -5, 10};
    Observable yvar{"yvar", -5, 10};

    // Fit parameter
    Variable mean1{"mean1", 1.5};
    Variable sigma1{"sigma1", 2.5};
    Variable mean2{"mean2", -1};
    Variable sigma2{"sigma2", 0.5};

    // GooPdf object
    GaussianPdf gausspdf1{"gausspdf1", xvar, mean1, sigma1};
    GaussianPdf gausspdf2{"gausspdf2", yvar, mean2, sigma2};

    ProdPdf product("product", {&gausspdf1, &gausspdf2});

    xvar = 1.2;
    yvar = -1.7;
    // std::cout << product.getValue(EvalFunc::Prob) << std::endl;
    // 0.085653187279569429
}
