#include <pybind11/pybind11.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/basic/BWPdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_BWPdf(py::module &m) {
    py::class_<BWPdf, GooPdf>(m, "BWPdf")
        .def(py::init<std::string, Variable*, Variable*, Variable*>())
        ;
}


