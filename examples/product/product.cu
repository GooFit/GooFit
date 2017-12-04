#include <goofit/Application.h>
#include <goofit/FitManager.h>
#include <goofit/PDFs/basic/ExpPdf.h>
#include <goofit/PDFs/combine/ProdPdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>

using namespace std;
using namespace GooFit;

int main(int argc, char **argv) {
    GooFit::Application app("Product example", argc, argv);

    GOOFIT_PARSE(app);

    Observable xvar{"xvar", 0, log(1 + RAND_MAX / 2)};
    Observable yvar{"yvar", 0, log(1 + RAND_MAX / 2)};
    vector<Observable> varList = {xvar, yvar};
    UnbinnedDataSet data{varList};

    for(int i = 0; i < 100000; ++i) {
        xvar.setValue(xvar.getUpperLimit() - log(1 + rand() / 2));
        yvar.setValue(yvar.getUpperLimit() - log(1 + rand() / 2));
        data.addEvent();
    }

    Variable alpha_x{"alpha_x", -2.4, 0.1, -10, 10};
    Variable alpha_y{"alpha_y", -1.1, 0.1, -10, 10};

    ExpPdf exp_x{"exp_x", xvar, alpha_x};
    ExpPdf exp_y{"exp_y", yvar, alpha_y};
    ProdPdf product{"product", {&exp_x, &exp_y}};

    product.setData(&data);
    FitManager fitter(&product);
    fitter.fit();

    return fitter;
}
