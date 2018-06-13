#include <goofit/Python.h>

#include <goofit/PDFs/physics/Amp3BodyBase.h>

using namespace GooFit;

void init_Amp3BodyBase(py::module &m) {
    py::class_<Amp3BodyBase, AmpNBodyBase>(m, "Amp3BodyBase")

        ;
}
