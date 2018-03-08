#include <pybind11/pybind11.h>

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <Minuit2/MinuitParameter.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ROOT::Minuit2;

void init_MinuitParameter(py::module &m) {
    py::class_<MinuitParameter>(m, "MinuitParameter")
        .def("Error", &MinuitParameter::Error)
        .def("Fix", &MinuitParameter::Fix)
        .def("GetName", &MinuitParameter::GetName)
        .def("HasLimits", &MinuitParameter::HasLimits)
        .def("HasLowerLimit", &MinuitParameter::HasLowerLimit)
        .def("HasUpperLimit", &MinuitParameter::HasUpperLimit)
        .def("IsConst", &MinuitParameter::IsConst)
        .def("IsFixed", &MinuitParameter::IsFixed)
        .def("LowerLimit", &MinuitParameter::LowerLimit)
        .def("Name", &MinuitParameter::Name)
        .def("Number", &MinuitParameter::Number)
        .def("Release", &MinuitParameter::Release)
        .def("RemoveLimits", &MinuitParameter::RemoveLimits)
        .def("SetError", &MinuitParameter::SetError)
        .def("SetLimits", &MinuitParameter::SetLimits)
        .def("SetLowerLimit", &MinuitParameter::SetLowerLimit)
        .def("SetName", &MinuitParameter::SetName)
        .def("SetUpperLimit", &MinuitParameter::SetUpperLimit)
        .def("SetValue", &MinuitParameter::SetValue)
        .def("UpperLimit", &MinuitParameter::UpperLimit)
        .def("Value", &MinuitParameter::Value);
}
