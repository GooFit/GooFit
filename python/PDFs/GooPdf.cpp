
#include <pybind11/pybind11.h>

#include <goofit/PDFs/GooPdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_GooPdf(py::module &m) {
    py::class_<GooPdf, PdfBase>(m, "GooPdf")
        ;

}

