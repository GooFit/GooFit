#include <pybind11/pybind11.h>

#include <goofit/PDFs/combine/ConvolutionPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_ConvolutionPdf(py::module &m) {
    py::class_<ConvolutionPdf, GooPdf>(m, "ConvolutionPdf")
        .def(py::init<std::string, Variable *, GooPdf *, GooPdf *>(),
             py::keep_alive<1,4>(),
             py::keep_alive<1,5>())
        .def(py::init<std::string, Variable *, GooPdf *, GooPdf *, unsigned int>(),
             py::keep_alive<1,4>(),
             py::keep_alive<1,5>());
}
