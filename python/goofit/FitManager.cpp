#include <goofit/Python.h>

#include <pybind11/iostream.h>
#include <pybind11/stl.h>

#include <Minuit2/MnUserParameters.h>
#include <Minuit2/FCNBase.h>

#include <goofit/PdfBase.h>
#include <goofit/PyProps.h>
#include <goofit/fitting/FitManagerMinuit2.h>
#include <goofit/fitting/Params.h>

#include <iostream>

using namespace GooFit;

void init_FitManager(py::module &m) {
    py::class_<FCN, Minuit2::FCNBase>(m, "FCN").def(py::init<Params &>(), "params"_a);

    py::class_<Params, Minuit2::MnUserParameters>(m, "Params")
        .def(py::init<PdfBase &>(), "pdf"_a)
        .def("SetGooFitParams", &Params::SetGooFitParams, "input"_a)
        .def("size", &Params::size)
        .def("make_minuit_vector", &Params::make_minuit_vector)

        // clang-format off
        ADD_PROP_WO(record, set_record, Params)
        ADD_PROP_RO(recorded, get_recorded, Params)
        // clang-format on
        ;

    py::class_<FitManagerMinuit2>(m, "FitManager")
        .def(py::init<PdfBase *>())
        .def("fit", &FitManagerMinuit2::fit, py::call_guard<py::scoped_ostream_redirect>())
        .def("getMinosErrors", &FitManagerMinuit2::getMinosErrors)
        .def("getMnScan", &FitManagerMinuit2::getMnScan)
        .def("__int__", &FitManagerMinuit2::operator int)
        .def("setMinos", &FitManagerMinuit2::setMinos)
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
