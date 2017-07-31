#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/physics/TruthResolution_Aux.h>

using namespace GooFit;
namespace py = pybind11;

void init_TruthResolution(py::module &m) {
    py::class_<TruthResolution, MixingTimeResolution>(m, "TruthResolution")
            .def(py::init<>())

            ;
}









