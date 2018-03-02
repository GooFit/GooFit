#pragma once

#include <goofit/PDFs/GooPdf.h>
#include <goofit/fitting/FitManagerMinuit1.h>

bool test_fitter(GooFit::GooPdf *pdf) {
    GooFit::FitManagerMinuit1 fitter{pdf};
    fitter.setVerbosity(2);
    fitter.fit();
    return bool(fitter);
}
