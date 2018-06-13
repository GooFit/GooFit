#include <goofit/Python.h>

#include <goofit/PDFs/physics/AmpNBodyBase.h>

using namespace GooFit;

void init_AmpNBodyBase(py::module &m) {
    py::class_<AmpNBodyBase, CombinePdf>(m, "AmpNBodyBase")

        ;
}
