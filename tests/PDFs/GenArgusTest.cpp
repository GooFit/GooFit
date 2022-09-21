#include "Pdf1D.h"
#include <goofit/Catch.h>

#include <goofit/PDFs/basic/ArgusPdf.h>

TEST_CASE_METHOD(PdfTest1D, "Argus generate and fit", "[gen][fit][argus]") {
    fptype answer;

    SECTION("Set to 1") { answer = 1.0; }
    SECTION("Set to 2") { answer = 2.0; }

    Variable c{"c", xvar.getUpperLimit()};
    Variable xi{"xi", answer, 0.1, 0.0, 10.0};

    set_pdf<ArgusPdf>(c, xi, true);

    // Generate toy events
    fill();

    // Change the values to make the fit interesting
    xi = 1.5;

    INFO(fit());
    REQUIRE(result());

    int n      = 10000;
    double err = xi.getValue() / sqrt(n);
    // TODO: Is this reasonable?
    CHECK(xi.getError() < 5. * err);
    CHECK(xi.getValue() == Approx(answer).epsilon(3. / sqrt(n)));
}
