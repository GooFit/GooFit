#include <pybind11/pybind11.h>

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <Minuit2/MnUserParameterState.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ROOT::Minuit2;

void init_MnUserParameterState(py::module &m) {
    py::class_<MnUserParameterState>(m, "MnUserParameterState")
        .def("Value", (double(MnUserParameterState::*)(unsigned int) const) & MnUserParameterState::Value, "n"_a)
        .def("Error", (double(MnUserParameterState::*)(unsigned int) const) & MnUserParameterState::Error, "n"_a);
}
