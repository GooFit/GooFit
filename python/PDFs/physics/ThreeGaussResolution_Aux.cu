#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/PDFs/physics/ThreeGaussResolution_Aux.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_ThreeGaussResolution(py::module &m) {
    py::class_<ThreeGaussResolution, MixingTimeResolution>(m, "ThreeGaussResolution")
        .def(py::init<Variable *, Variable *, Variable *, Variable *, Variable *, Variable *, Variable *, Variable *>(),
             py::keep_alive<1, 2>(),
             py::keep_alive<1, 3>(),
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>(),
             py::keep_alive<1, 6>(),
             py::keep_alive<1, 7>(),
             py::keep_alive<1, 8>(),
             py::keep_alive<1, 9>());
}
