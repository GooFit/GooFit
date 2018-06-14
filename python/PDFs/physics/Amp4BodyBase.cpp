#include <goofit/Python.h>

#include <goofit/PDFs/physics/Amp4BodyBase.h>

using namespace GooFit;

void init_Amp4BodyBase(py::module &m) {
    py::class_<Amp4BodyBase, AmpNBodyBase>(m, "Amp4BodyBase")

        ;
}
