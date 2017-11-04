#include <pybind11/pybind11.h>

#include <goofit/PDFs/basic/ExpPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_ExpPdf(py::module &m) {
    py::class_<ExpPdf, GooPdf>(m, "ExpPdf")
        .def(py::init<std::string, Observable, Variable>(), py::keep_alive<1, 3>(), py::keep_alive<1, 4>())
        .def(py::init<std::string, Observable, Variable, Variable>(),
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>(),
             py::keep_alive<1, 3>());
}
