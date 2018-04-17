#include <pybind11/pybind11.h>

#include <goofit/PDFs/basic/ExpPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_ExpPdf(py::module &m) {
    py::class_<ExpPdf, GooPdf>(m, "ExpPdf")
        .def(py::init<std::string, Observable, Variable>(), "n", "_x", "alpha")
        .def(py::init<std::string, Observable, Variable, Variable>(), "n", "_x", "alpha", "offset");
}
