#include <pybind11/pybind11.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/combine/CompositePdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_CompositePdf(py::module &m) {
    py::class_<CompositePdf, GooPdf>(m, "CompositePdf")
        .def(py::init<std::string, PdfBase *, PdfBase *>())
        ;
}




