#include <pybind11/pybind11.h>

#include <pybind11/stl.h>

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnApplication.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ROOT::Minuit2;

void init_MnApplication(py::module &m) {
    py::class_<MnApplication>(m, "MnApplication")
        .def("__call__",
             &MnApplication::operator(),
             "Minimize the function, returns a function minimum",
             "maxfcn"_a    = 0,
             "tolerance"_a = 0.1);
}
