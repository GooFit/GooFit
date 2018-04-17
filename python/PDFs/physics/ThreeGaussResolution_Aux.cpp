#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/PDFs/physics/ThreeGaussResolution_Aux.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_ThreeGaussResolution(py::module &m) {
    py::class_<ThreeGaussResolution, MixingTimeResolution>(m, "ThreeGaussResolution")
        .def(py::init<Variable, Variable, Variable, Variable, Variable, Variable, Variable, Variable>(),
             "cf",
             "tf",
             "cb",
             "cs",
             "tb",
             "ts",
             "ob",
             "os");
}
