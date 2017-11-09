#include <pybind11/pybind11.h>

#include <goofit/PDFs/combine/ConvolutionPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_ConvolutionPdf(py::module &m) {
    py::class_<ConvolutionPdf, GooPdf>(m, "ConvolutionPdf")
        .def(py::init<std::string, Observable, GooPdf *, GooPdf *>(),
             "n",
             "_x",
             "model",
             "resolution",
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>())
        .def(py::init<std::string, Observable, GooPdf *, GooPdf *, unsigned int>(),
             "n",
             "_x",
             "model",
             "resolution",
             "numOthers",
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>());
}
