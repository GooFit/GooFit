#include <pybind11/pybind11.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/basic/JohnsonSUPdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_JohnsonSUPdf(py::module &m) {
    py::class_<JohnsonSUPdf, GooPdf>(m, "JohnsonSUPdf")
        .def(py::init<std::string, Variable*, Variable*, Variable*, Variable*, Variable*>())
        ;
}



