/*
04/05/2016 Christoph Hasse
DISCLAIMER:

This code is not sufficiently tested yet and still under heavy development!
See *.cu file for more details
*/

#pragma once

#include <goofit/PDFs/physics/Amplitude.h>
#include <goofit/PDFs/physics/lineshapes/Bugg.h>
#include <goofit/PDFs/physics/lineshapes/Bugg3.h>
#include <goofit/PDFs/physics/lineshapes/FOCUS.h>
#include <goofit/PDFs/physics/lineshapes/Flatte.h>
#include <goofit/PDFs/physics/lineshapes/GLASS.h>
#include <goofit/PDFs/physics/lineshapes/GSpline.h>
#include <goofit/PDFs/physics/lineshapes/LASS.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/lineshapes/NonRes.h>
#include <goofit/PDFs/physics/lineshapes/One.h>
#include <goofit/PDFs/physics/lineshapes/RBW.h>
#include <goofit/PDFs/physics/lineshapes/SBW.h>
#if GOOFIT_KMATRIX
#include <goofit/PDFs/physics/lineshapes/kMatrix.h>
#endif
