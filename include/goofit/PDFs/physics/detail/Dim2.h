/*
04/18/2016 Christoph Hasse

This file includes a converter from 16 values (4 4-momenta) to (the common set of) 5 parameters.
See UserUtils for a handy way to convert for tests.

*/

#pragma once

#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <mcbooster/GContainers.h>
#include <mcbooster/GFunctional.h>
#include <mcbooster/GTypes.h>
#include <mcbooster/Vector4R.h>

namespace GooFit {

struct Dim2 : public mcbooster::IFunctionArray {
    Dim2() { dim = 2; }

    __host__ __device__ void
    operator()(const mcbooster::GInt_t n, mcbooster::Vector4R **particles, mcbooster::GReal_t *variables) override {
        mcbooster::Vector4R p1 = *particles[0];
        mcbooster::Vector4R p2 = *particles[1];
        mcbooster::Vector4R p3 = *particles[2];

        mcbooster::Vector4R p12 = p1 + p2;
        mcbooster::Vector4R p23 = p2 + p3;
        mcbooster::Vector4R p13 = p1 + p3;

        variables[0] = p12.mass2();
        variables[1] = p23.mass2();
        variables[2] = p13.mass2();
    }
};

} // namespace GooFit
