#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/combine/ProdPdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_ProdPdf(py::module &m) {
    py::class_<ProdPdf, GooPdf>(m, "ProdPdf")
        .def(py::init<std::string, std::vector<PdfBase *>>())
        ;
}







