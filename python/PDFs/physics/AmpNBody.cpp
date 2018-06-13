#include <goofit/Python.h>

#include <goofit/PDFs/physics/AmpNBody.h>

using namespace GooFit;

void init_AmpNBody(py::module &m) {
    py::class_<AmpNBody, CombinePdf>(m, "AmpNBody")

        ;
}
