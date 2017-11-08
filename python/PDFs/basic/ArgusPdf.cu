#include <pybind11/pybind11.h>

#include <goofit/PDFs/basic/ArgusPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_ArgusPdf(py::module &m) {
    py::class_<ArgusPdf, GooPdf>(m, "ArgusPdf")
        .def(py::init<std::string, Observable, Variable, Variable, bool>(),
             "n",
             "_x",
             "m",
             "s",
             "upper",
             py::keep_alive<1, 3>(),
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>())
        .def(py::init<std::string, Observable, Variable, Variable, bool, Variable>(),
             "n",
             "_x",
             "m",
             "s",
             "upper",
             "power",
             py::keep_alive<1, 3>(),
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>());

}
