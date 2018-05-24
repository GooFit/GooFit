#include "Pdf1D.h"
#include <goofit/Catch.h>

#include <goofit/PDFs/basic/GaussianPdf.h>

TEST_CASE_METHOD(PdfTest1D, "Gaussian", "[fit][gaussian]") {
    xvar.setLowerLimit(-10);
    xvar.setUpperLimit(10);

    fptype mean_ans, sigma_ans;

    SECTION("Set to 1, 2") {
        mean_ans  = 1.0;
        sigma_ans = 2.0;
    }
    SECTION("Set to 2, 1") {
        mean_ans  = 2.0;
        sigma_ans = 1.0;
    }

    Variable mean{"mean", mean_ans, -10, 10};
    Variable sigma{"sigma", sigma_ans, 0.0, 10.0};

    set_pdf<GaussianPdf>(mean, sigma);

    // Generate toy events
    fill();

    // Change the values to make the fit interesting
    mean  = 0.0;
    sigma = 1.5;

    INFO(fit());
    REQUIRE(result());

    CHECK(mean.getError() < .01);
    CHECK(mean.getValue() == Approx(mean_ans).epsilon(.01));

    CHECK(sigma.getError() < .01);
    CHECK(sigma.getValue() == Approx(sigma_ans).epsilon(.01));
}
