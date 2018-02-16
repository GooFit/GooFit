#include <iostream>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>

#include <goofit/FitManager.h>
#include <goofit/PdfBase.h>

#include "props.h"

namespace py = pybind11;
using namespace GooFit;

void init_FitManager(py::module &m) {
    py::class_<FitManager>(m, "FitManager")
        .def(py::init<PdfBase *>())
        .def("fit", &FitManager::fit, py::call_guard<py::scoped_ostream_redirect>())
        .def("__int__", &FitManager::operator int)
        .def("__bool__", &FitManager::operator bool)    // Py3
        .def("__nonzero__", &FitManager::operator bool) // Py2
        // clang-format off
        ADD_PROP(verbosity, getVerbosity, setVerbosity, FitManager)
        ADD_PROP_WO(max_calls, setMaxCalls, FitManager)
        // clang-format on
        ;
}
