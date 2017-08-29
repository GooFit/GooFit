#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

#include <Minuit2/FunctionMinimum.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ROOT::Minuit2;


void init_FunctionMinimum(py::module& m) {

    py::class_<MnUserParameterState>(m, "MnUserParameterState")
        .def("Value", (double (MnUserParameterState::*)(unsigned int) const) &MnUserParameterState::Value, "n"_a)
        .def("Error", (double (MnUserParameterState::*)(unsigned int) const) &MnUserParameterState::Error, "n"_a)
    ;

    py::class_<MnUserParameters>(m, "MnUserParameters")
    ;


    py::class_<MnUserCovariance>(m, "MnUserCovariance")
    ;

    py::class_<FunctionMinimum>(m, "FunctionMinimum")
        .def("UserState", &FunctionMinimum::UserState)
        .def("UserParameters", &FunctionMinimum::UserParameters)
        .def("UserCovariance", &FunctionMinimum::UserCovariance)
    ;
}

