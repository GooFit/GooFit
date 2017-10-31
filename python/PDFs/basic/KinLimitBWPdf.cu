#include <pybind11/pybind11.h>

#include <goofit/PDFs/basic/KinLimitBWPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_KinLimitBWPdf(py::module &m) {
    py::class_<KinLimitBWPdf, GooPdf>(m, "KinLimitBWPdf")
        .def(py::init<std::string, Variable *, Variable *, Variable *>(),
             py::keep_alive<1, 3>(),
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>());
}
