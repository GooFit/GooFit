#pragma once

#include <goofit/fitting/FitManagerMinuit1.h>
#include <goofit/PDFs/GooPdf.h>

bool test_fitter(GooFit::GooPdf* pdf) {
    GooFit::FitManagerMinuit1 fitter{pdf};
    fitter.setVerbosity(2);
    fitter.fit();
    return bool(fitter);
}

