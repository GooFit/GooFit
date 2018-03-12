#pragma once

#include <goofit/Version.h>

#if GOOFIT_ROOT_FOUND

#include <RVersion.h>
#include <TCanvas.h>
#include <TStyle.h>

#endif

namespace GooFit {

#if GOOFIT_ROOT_FOUND

/// A default set of ROOT styles for pretty GooFit plots
void setROOTStyle() {
    gStyle->SetNumberContours(512);
    gStyle->SetOptStat(0);
    gStyle->SetPadRightMargin(.13);

    if(ROOT_VERSION_CODE < ROOT_VERSION(6, 6, 0))
        gStyle->SetPalette(kRainBow, 0);
    else
        gStyle->SetPalette(kViridis, 0);
}

#endif

} // namespace GooFit
