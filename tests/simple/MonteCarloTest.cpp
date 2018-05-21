#include <goofit/Catch.h>

#include <goofit/FitManager.h>
#include <goofit/PDFs/basic/ExpPdf.h>
#include <goofit/PDFs/basic/GaussianPdf.h>
#include <goofit/PDFs/combine/ProdPdf.h>
#include <goofit/UnbinnedDataSet.h>

#include <goofit/Variable.h>

using namespace GooFit;

TEST_CASE("Monte Carlo Gaussian generate and fit", "[simple][gauss][fit][gen]") {
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
    pdf.fillMCDataSimple(1000, /* seed */ 27);

    fmt::print("Data set entries: {} \n", data.getNumEvents());

    // Change the fitting parameters
    mu.setValue(1.2);
    sigma.setValue(.4);

    FitManager fitter{&pdf};
    fitter.fit();

    CHECK(fitter);
    CHECK(mu.getError() < .1);
    CHECK(sigma.getError() < .1);
    CHECK(mu.getValue() == Approx(.5).margin(.3));
    CHECK(sigma.getValue() == Approx(.8).margin(.3));
}
