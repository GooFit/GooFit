#include <pybind11/pybind11.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/combine/ConvolutionPdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_ConvolutionPdf(py::module &m) {
    py::class_<ConvolutionPdf, GooPdf>(m, "ConvolutionPdf")
        .def(py::init<std::string, Variable *, GooPdf *, GooPdf *>())
        .def(py::init<std::string, Variable *, GooPdf *, GooPdf *, unsigned int>())
        ;
}




