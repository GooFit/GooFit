#include <pybind11/pybind11.h>

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <Minuit2/MnUserParameters.h>

#include <string>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ROOT::Minuit2;

void init_MnUserParameters(py::module &m) {
    py::class_<MnUserParameters>(m, "MnUserParameters")
        .def(py::init<>())
        .def("Add", (bool (MnUserParameters::*)(const std::string &, double)) & MnUserParameters::Add)
        .def("Add", (bool (MnUserParameters::*)(const std::string &, double, double)) & MnUserParameters::Add)
        .def("Add",
             (bool (MnUserParameters::*)(const std::string &, double, double, double, double)) & MnUserParameters::Add)
        .def("Parameter", &MnUserParameters::Parameter, "Get a parameter by number")
        .def("Parameters", &MnUserParameters::Parameters, "Get all the parameters");
}
