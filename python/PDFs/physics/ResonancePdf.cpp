#include <pybind11/pybind11.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/physics/ResonancePdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_ResonancePdf(py::module &m) {
    py::class_<ResonancePdf, GooPdf>(m, "ResonancePdf")
        .def(py::init<std::string, Variable *, Variable *, Variable *, Variable *, unsigned int, unsigned int>())
        .def(py::init<std::string, Variable *, Variable *, unsigned int, Variable *, Variable *, unsigned int>())
        .def(py::init<std::string, Variable *, Variable *, Variable *, unsigned int, Variable *, unsigned int>())
        .def(py::init<std::string, Variable *, Variable *>())
        .def(py::init<std::string, Variable *, Variable *, Variable *, Variable *, unsigned int>())

        ;
}








