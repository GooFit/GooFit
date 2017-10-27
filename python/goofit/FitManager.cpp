#include <iostream>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>

#include <goofit/FitManager.h>
#include <goofit/PdfBase.h>

namespace py = pybind11;
using namespace GooFit;

void init_FitManager(py::module &m) {
    py::class_<FitManager>(m, "FitManager")
        .def(py::init<PdfBase *>())
        .def("fit", &FitManager::fit, py::call_guard<py::scoped_ostream_redirect>())
        .def("__int__", &FitManager::operator int)
        ;
}
