#include <goofit/Catch.h>

#include <goofit/DataSet.h>
#include <goofit/PDFs/basic/GaussianPdf.h>
#include <goofit/Workspace.h>

using namespace GooFit;
using Catch::Matchers::Contains;

TEST_CASE("Checking Workspace printout", "[simple][workspace]") {
    Workspace w;

    w.add_obs(Observable{"xvar", -10.0, 10.0});

    w.add_var(Variable{"mu", 1.0});
    w.add_var(Variable{"sigma", 2.0});

    w.add_pdf(new GaussianPdf{"gauss", *w.obs("xvar"), *w.var("mu"), *w.var("sigma")});

    w.add_dataset(new UnbinnedDataSet{*w.obs("xvar")});

    std::cout << w.Print();
}
