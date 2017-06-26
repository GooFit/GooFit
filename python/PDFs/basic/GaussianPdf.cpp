#include <pybind11/pybind11.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/basic/GaussianPdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_GaussianPdf(py::module &m) {
    py::class_<GaussianPdf, GooPdf>(m, "GaussianPdf")
        .def(py::init<std::string, Variable*, Variable*, Variable*>())
        ;
}


