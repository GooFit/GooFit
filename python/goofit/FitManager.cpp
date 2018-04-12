#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/PdfBase.h>
#include <goofit/fitting/FitManagerMinuit2.h>

#include <iostream>

#include "props.h"

namespace py = pybind11;
using namespace GooFit;

void init_FitManager(py::module &m) {
    py::class_<FCN>(m, "FCN");

    py::class_<Params>(m, "Params")
        // clang-format off
        ADD_PROP_WO(record, set_record, Params)
        ADD_PROP_RO(recorded, get_recorded, Params)
        // clang-format on
        ;

    py::class_<FitManagerMinuit2>(m, "FitManager")
        .def(py::init<PdfBase *>())
        .def("fit", &FitManagerMinuit2::fit, py::call_guard<py::scoped_ostream_redirect>())
        .def("__int__", &FitManagerMinuit2::operator int)
        .def("__bool__", &FitManagerMinuit2::operator bool)    // Py3
        .def("__nonzero__", &FitManagerMinuit2::operator bool) // Py2

        // clang-format off
        ADD_PROP(verbosity, getVerbosity, setVerbosity, FitManagerMinuit2)
        ADD_PROP_WO(max_calls, setMaxCalls, FitManagerMinuit2)

        ADD_PROP_RO(params, getParams, FitManagerMinuit2)
        ADD_PROP_RO(fcn, getFCN, FitManagerMinuit2)
        // clang-format on
        ;
}
