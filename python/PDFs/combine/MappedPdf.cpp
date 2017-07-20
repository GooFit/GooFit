#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/combine/MappedPdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_MappedPdf(py::module &m) {
    py::class_<MappedPdf, GooPdf>(m, "MappedPdf")
        .def(py::init<std::string, GooPdf *,std::vector<GooPdf *>>())
        ;
}






