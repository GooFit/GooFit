#include <pybind11/pybind11.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/basic/ArgusPdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_ArgusPdf(py::module &m) {
    py::class_<ArgusPdf, GooPdf>(m, "ArgusPdf")
        .def(py::init<std::string, Variable*, Variable*, Variable*, bool>())
        .def(py::init<std::string, Variable*, Variable*, Variable*, bool, Variable*>())
        ;
}

