#include <pybind11/pybind11.h>

#include <goofit/PDFs/basic/ScaledGaussianPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_ScaledGaussianPdf(py::module &m) {
    py::class_<ScaledGaussianPdf, GooPdf>(m, "ScaledGaussianPdf")
        .def(py::init<std::string, Variable *, Variable *, Variable *, Variable *, Variable *>(),
             py::keep_alive<1,3>(),
             py::keep_alive<1,4>(),
             py::keep_alive<1,5>(),
             py::keep_alive<1,6>(),
             py::keep_alive<1,7>());
}
