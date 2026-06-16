#include <goofit/Catch.h>

#include <goofit/PDFs/basic/GaussianPdf.h>
#include <goofit/PDFs/combine/AddPdf.h>
#include <goofit/UnbinnedDataSet.h>

#include <goofit/Variable.h>

using namespace GooFit;

// Regression test: a non-extended AddPdf f*A + (1-f)*B must normalize *every*
// component, including the last one. Picking two Gaussians with different widths
// makes their normalization factors differ, so a missing factor on the last term
// is detectable. The combined probability must equal the weighted sum of the
// (already normalized) component probabilities at every point.
TEST_CASE("AddPdf component normalization", "[basic]") {
    Observable xvar{"xvar", -5, 10};

    Variable mean{"mean", 1.5};
    Variable sigmaA{"sigmaA", 2.5};
    Variable sigmaB{"sigmaB", 0.8};

    GaussianPdf gaussA{"gaussA", xvar, mean, sigmaA};
    GaussianPdf gaussB{"gaussB", xvar, mean, sigmaB};

    Variable frac{"frac", 0.6};

    // Non-extended two-component constructor: frac * gaussA + (1 - frac) * gaussB.
    AddPdf sum{"sum", frac, &gaussA, &gaussB};

    auto dataset = sum.makeGrid();
    sum.setData(&dataset);

    auto vv = sum.getCompProbsAtDataPoints();

    REQUIRE(vv.size() == 3); // combined, gaussA, gaussB
    REQUIRE(vv[0].size() == vv[1].size());
    REQUIRE(vv[0].size() == vv[2].size());

    const fptype f = 0.6;
    for(size_t i = 0; i < vv[0].size(); i += 7) {
        CHECK(vv[0][i] == Approx(f * vv[1][i] + (1 - f) * vv[2][i]));
    }
}
