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
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasColor(10);
    gStyle->SetFrameFillColor(10);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetPadColor(0);
    gStyle->SetTitleColor(1);
    gStyle->SetStatColor(0);
    gStyle->SetFillColor(0);
    gStyle->SetFuncWidth(1);
    gStyle->SetLineWidth(1);
    gStyle->SetLineColor(1);
    gStyle->SetNumberContours(256);
    gStyle->SetOptStat(0);

    if(ROOT_VERSION_CODE < ROOT_VERSION(6, 6, 0))
        gStyle->SetPalette(kRainBow, 0);
    else
        gStyle->SetPalette(kViridis, 0);
}

#endif

} // namespace GooFit
