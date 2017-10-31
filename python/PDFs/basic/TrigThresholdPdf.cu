#include <pybind11/pybind11.h>

#include <goofit/PDFs/basic/TrigThresholdPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_TrigThresholdPdf(py::module &m) {
    py::class_<TrigThresholdPdf, GooPdf>(m, "TrigThresholdPdf")
        .def(py::init<std::string, Variable *, Variable *, Variable *, Variable *>(),
             py::keep_alive<1, 3>(),
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>(),
             py::keep_alive<1, 6>())
        .def(py::init<std::string, Variable *, Variable *, Variable *, Variable *, bool>(),
             py::keep_alive<1, 3>(),
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>(),
             py::keep_alive<1, 6>())
        .def(py::init<std::string, Variable *, Variable *, Variable *, Variable *, Variable *, Variable *, bool>(),
             py::keep_alive<1, 3>(),
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>(),
             py::keep_alive<1, 6>(),
             py::keep_alive<1, 7>(),
             py::keep_alive<1, 8>());
}
