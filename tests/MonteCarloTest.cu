#include <gtest/gtest.h>

#include <goofit/FitManager.h>
#include <goofit/PDFs/basic/ExpPdf.h>
#include <goofit/PDFs/basic/GaussianPdf.h>
#include <goofit/PDFs/combine/ProdPdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/detail/MonteCarlo.h>

#include <goofit/Variable.h>

using namespace GooFit;

TEST(MonteCarlo, GaussianGenAndFit) {
    // Independent variable.
    Observable var{"var", -10, 10};
    var.setNumBins(100);

    // Data set
    UnbinnedDataSet data(var);

    // Fit parameter
    Variable mu{"mu", .5, 0.1, -10, 10};
    Variable sigma{"sigma", .8, 0.1, -10, 10};

    // GooPdf object
    GaussianPdf pdf{"pdf", var, mu, sigma};
    pdf.setData(&data);

    // Generate toy events.
    fillDataSetMC1D(pdf, var, 1000, /* seed */ 27);

    fmt::print("Data set entries: {} \n", data.getNumEvents());

    // Change the fitting parameters
    mu.setValue(1.2);
    sigma.setValue(.4);

    FitManager fitter{&pdf};
    fitter.fit();

    EXPECT_TRUE(fitter);
    EXPECT_LT(mu.getError(), .1);
    EXPECT_LT(sigma.getError(), .1);
    EXPECT_NEAR(.5, mu.getValue(), mu.getError() * 3);
    EXPECT_NEAR(.8, sigma.getValue(), sigma.getError() * 3);
}
