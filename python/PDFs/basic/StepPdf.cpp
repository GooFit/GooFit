#include <pybind11/pybind11.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/basic/StepPdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_StepPdf(py::module &m) {
    py::class_<StepPdf, GooPdf>(m, "StepPdf")
        .def(py::init<std::string, Variable *, Variable *>())
        ;
}









