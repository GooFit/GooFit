#include <pybind11/pybind11.h>

#include <goofit/PDFs/basic/ArgusPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_ArgusPdf(py::module &m) {
    py::class_<ArgusPdf, GooPdf>(m, "ArgusPdf")
        .def(py::init<std::string, Variable *, Variable *, Variable *, bool>(),
             py::keep_alive<1, 3>(),
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>())
        .def(py::init<std::string, Variable *, Variable *, Variable *, bool, Variable *>(),
             py::keep_alive<1, 3>(),
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>());
}
