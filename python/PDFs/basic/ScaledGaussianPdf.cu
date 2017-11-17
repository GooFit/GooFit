#include <pybind11/pybind11.h>

#include <goofit/PDFs/basic/ScaledGaussianPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_ScaledGaussianPdf(py::module &m) {
    py::class_<ScaledGaussianPdf, GooPdf>(m, "ScaledGaussianPdf")
        .def(
            py::init<std::string, Observable, Variable, Variable, Variable, Variable>(), "n", "_x", "m", "s", "d", "e");
}
