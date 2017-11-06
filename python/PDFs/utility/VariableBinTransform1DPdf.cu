#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/PDFs/utility/VariableBinTransform1DPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_VariableBinTransform1DPdf(py::module &m) {
    py::class_<VariableBinTransform1DPdf, GooPdf>(m, "VariableBinTransform1DPdf")
        .def(py::init<std::string, Observable, std::vector<fptype>>());
}
