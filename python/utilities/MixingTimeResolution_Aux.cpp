
#include <pybind11/pybind11.h>

#include <goofit/PDFs/physics/MixingTimeResolution_Aux.h>

using namespace GooFit;
namespace py = pybind11;

void init_MixingTimeResolution(py::module &m) {
    py::class_<MixingTimeResolution>(m, "MixingTimeResolution")
        .def("initIndex", &MixingTimeResolution::initIndex)
    ;
}

