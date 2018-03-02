#include <gtest/gtest.h>

#include <goofit/FitManager.h>
#include <goofit/PDFs/basic/ExpPdf.h>
#include <goofit/PDFs/basic/GaussianPdf.h>
#include <goofit/PDFs/combine/ProdPdf.h>
#include <goofit/UnbinnedDataSet.h>

#include <goofit/Variable.h>

#include <iostream>
#include <sys/time.h>
#include <sys/times.h>

#include <random>

using namespace GooFit;

TEST(PDFComps, KnownNormalize) {
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

    EXPECT_FLOAT_EQ(exppdf.integrate(0, 10), 0.5);
    EXPECT_FLOAT_EQ(exppdf.normalize(), 0.5);
    EXPECT_FLOAT_EQ(exppdf.normalize(), 0.5); // Just verifying that it does not crash

    FitManager fitter{&exppdf};
    fitter.setVerbosity(0);
    fitter.fit();

    EXPECT_FLOAT_EQ(exppdf.normalize(), 0.665099);
    EXPECT_FLOAT_EQ(exppdf.normalize(), 0.665099);
    ; // Just verifying that it does not crash
    EXPECT_FLOAT_EQ(exppdf.normalize(), exppdf.integrate(0, 10));
}

/*
 import numpy as np
 from scipy.stats import norm
 bin_centers = np.linspace(-5,10,101) + .15/2
 vv = norm.pdf(bin_centers,1.5,2.5)
 for i in range(0,80,20):
     print(i, bin_centers[i], vv[i])
 */

TEST(PDFComps, OneDGrid) {
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

    EXPECT_FLOAT_EQ(dataset.getValue(xvar, 0), -4.925);
    EXPECT_FLOAT_EQ(dataset.getValue(xvar, 20), -1.925);
    EXPECT_FLOAT_EQ(dataset.getValue(xvar, 40), 1.075);
    EXPECT_FLOAT_EQ(dataset.getValue(xvar, 60), 4.075);

    // Compute probabilities at points
    auto vv = gausspdf.getCompProbsAtDataPoints();

    // Check
    ASSERT_EQ(vv.size(), 1);
    ASSERT_EQ(vv[0].size(), 100);

    EXPECT_FLOAT_EQ(vv[0][0], 0.00587129964482);
    EXPECT_FLOAT_EQ(vv[0][20], 0.0624318782242);
    EXPECT_FLOAT_EQ(vv[0][40], 0.157287605852);
    EXPECT_FLOAT_EQ(vv[0][60], 0.0938855055588);
}

/*
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

 plt.contourf(X,Y,vv,40);
 */

TEST(PDFComps, TwoDGrid) {
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

    EXPECT_FLOAT_EQ(dataset.getValue(xvar, 0), -4.925);
    EXPECT_FLOAT_EQ(dataset.getValue(xvar, 20), -1.925);
    EXPECT_FLOAT_EQ(dataset.getValue(xvar, 40), 1.075);
    EXPECT_FLOAT_EQ(dataset.getValue(xvar, 60), 4.075);

    EXPECT_FLOAT_EQ(dataset.getValue(yvar, 100 * 0), -4.925);
    EXPECT_FLOAT_EQ(dataset.getValue(yvar, 100 * 20), -1.925);
    EXPECT_FLOAT_EQ(dataset.getValue(yvar, 100 * 40), 1.075);
    EXPECT_FLOAT_EQ(dataset.getValue(yvar, 100 * 60), 4.075);

    // Compute probabilities at points
    auto vv = product.getCompProbsAtDataPoints();

    // Check
    ASSERT_EQ(vv.size(), 3);
    ASSERT_EQ(vv[0].size(), 10000);
    ASSERT_EQ(vv[1].size(), 10000);
    ASSERT_EQ(vv[2].size(), 10000);

    EXPECT_FLOAT_EQ(vv[1][0], 0.00587129964482);
    EXPECT_FLOAT_EQ(vv[1][20], 0.0624318782242);
    EXPECT_FLOAT_EQ(vv[1][40], 0.157287605852);
    EXPECT_FLOAT_EQ(vv[1][60], 0.0938855055588);

    // EXPECT_FLOAT_EQ(vv[0][0], 7.53698529533e-12);
    // EXPECT_FLOAT_EQ(vv[0][20*100+20], 0.00579245632235);
    // EXPECT_FLOAT_EQ(vv[0][20*100+60], 0.0160636095436);
}

TEST(PDFComps, OneDEval) {
    // Independent variable.
    Observable xvar{"xvar", -5, 10};

    // Fit parameter
    Variable mean{"mean", 1.5};
    Variable sigma{"sigma", 2.5};

    // GooPdf object
    GaussianPdf gauss{"gausspdf", xvar, mean, sigma};

    auto v = gauss.evaluateAtPoints(xvar);

    EXPECT_FLOAT_EQ(v[0], 0.00587129964482);
    EXPECT_FLOAT_EQ(v[20], 0.0624318782242);
    EXPECT_FLOAT_EQ(v[40], 0.157287605852);
    EXPECT_FLOAT_EQ(v[60], 0.0938855055588);

    EXPECT_FLOAT_EQ(gauss.normalize(), 6.26657);

    xvar = 1.2;
    EXPECT_FLOAT_EQ(gauss.getValue(EvalFunc::Prob), 0.15843208471746245);
    // The default is to ignore the normalization factors, so I'm passing EvalFunc::Prob

    xvar = 7.2;
    EXPECT_FLOAT_EQ(gauss.getValue(EvalFunc::Prob), 0.0118618339389365);
}

TEST(PDFComps, TwoDEval) {
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
