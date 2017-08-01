#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/basic/BinTransformPdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_BinTransformPdf(py::module &m) {
    py::class_<BinTransformPdf, GooPdf>(m, "BinTransformPdf")
        .def(py::init<std::string, std::vector<Variable *>, std::vector<fptype>, std::vector<fptype>, std::vector<int>>())
        ;
}


