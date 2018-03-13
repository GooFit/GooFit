#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <sstream>

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnPrint.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ROOT::Minuit2;

void init_FunctionMinimum(py::module &m) {
    py::class_<FunctionMinimum>(m, "FunctionMinimum")
        .def("UserState", &FunctionMinimum::UserState)
        .def("UserParameters", &FunctionMinimum::UserParameters)
        .def("UserCovariance", &FunctionMinimum::UserCovariance)

        .def("Fval", &FunctionMinimum::Fval)
        .def("Edm", &FunctionMinimum::Edm)
        .def("NFcn", &FunctionMinimum::NFcn, "Number of function calls")

        .def("Up", &FunctionMinimum::Up)
        .def("IsValid", &FunctionMinimum::IsValid)
        .def("HasValidParameters", &FunctionMinimum::HasValidParameters)
        .def("HasValidCovariance", &FunctionMinimum::HasValidCovariance)
        .def("HasAccurateCovar", &FunctionMinimum::HasAccurateCovar)
        .def("HasPosDefCovar", &FunctionMinimum::HasPosDefCovar)
        .def("HasMadePosDefCovar", &FunctionMinimum::HasMadePosDefCovar)
        .def("HesseFailed", &FunctionMinimum::HesseFailed)
        .def("HasCovariance", &FunctionMinimum::HasCovariance)
        .def("IsAboveMaxEdm", &FunctionMinimum::IsAboveMaxEdm)
        .def("HasReachedCallLimit", &FunctionMinimum::HasReachedCallLimit)

        .def("__str__", [](const FunctionMinimum &self) {
            std::stringstream os;
            os << self;
            return os.str();
        });
}
