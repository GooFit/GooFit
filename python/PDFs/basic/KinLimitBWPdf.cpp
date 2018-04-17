#include <pybind11/pybind11.h>

#include <goofit/PDFs/basic/KinLimitBWPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_KinLimitBWPdf(py::module &m) {
    py::class_<KinLimitBWPdf, GooPdf>(m, "KinLimitBWPdf")
        .def(py::init<std::string, Observable, Variable, Variable>(), "n", "_x", "m", "s");
}
