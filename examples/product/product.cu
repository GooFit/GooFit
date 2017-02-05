#include "goofit/Application.h"
#include "goofit/Variable.h"
#include "goofit/FitManager.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/PDFs/ExpPdf.h"
#include "goofit/PDFs/ProdPdf.h"

using namespace std;

int main(int argc, char** argv) {

    GooFit::Application app("Product example", argc, argv);

    try {
        app.run();
    } catch (const GooFit::Error &e) {
        return app.exit(e);
    }

    Variable xvar{"xvar", 0, log(1+RAND_MAX/2)};
    Variable yvar{"yvar", 0, log(1+RAND_MAX/2)};
    vector<Variable*> varList = {&xvar, &yvar};
    UnbinnedDataSet data{varList};

    for(int i = 0; i < 100000; ++i) {
        xvar.value = xvar.upperlimit - log(1+rand()/2);
        yvar.value = yvar.upperlimit - log(1+rand()/2);
        data.addEvent();
    }

    Variable alpha_x{"alpha_x", -2.4, 0.1, -10, 10};
    Variable alpha_y{"alpha_y", -1.1, 0.1, -10, 10};

    ExpPdf exp_x{"exp_x", &xvar, &alpha_x};
    ExpPdf exp_y{"exp_y", &yvar, &alpha_y};
    ProdPdf product{"product", {&exp_x, &exp_y}};
    
    product.setData(&data);
    FitManager fitter(&product);
    fitter.fit();

    return 0;
}
