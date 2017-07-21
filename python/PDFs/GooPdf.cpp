
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/PDFs/GooPdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_GooPdf(py::module &m) {
    py::class_<GooPdf, PdfBase>(m, "GooPdf")
        .def("makeGrid", &GooPdf::makeGrid)
        .def("evaluateAtPoints", &GooPdf::evaluateAtPoints)
    ;
}
