#pragma once

#include <goofit/PDFs/GooPdf.h>
#include <goofit/fitting/FitManagerMinuit2.h>

#include <goofit/Version.h>

#if GOOFIT_ROOT_FOUND
#include <goofit/fitting/FitManagerMinuit1.h>

bool test_fitter_minuit1(GooFit::GooPdf *pdf) {
    GooFit::FitManagerMinuit1 fitter{pdf};
    fitter.setVerbosity(2);
    fitter.fit();
    return bool(fitter);
}

#endif

bool test_fitter(GooFit::GooPdf *pdf) {
    GooFit::FitManagerMinuit2 fitter{pdf};
    fitter.setVerbosity(2);
    fitter.fit();
    return bool(fitter);
}
