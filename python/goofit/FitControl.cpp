#include <iostream>

#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/FitControl.h>
#include <Minuit2/FunctionMinimum.h>

namespace py = pybind11;
using namespace GooFit;
using namespace pybind11::literals;

void init_FitControl(py::module &m) {
    py::class_<FitControl>(m, "FitControl")
             .def("binnedFit", &FitControl::binnedFit)
             .def("binErrors", &FitControl::binErrors)
             .def("metricIsPdf", &FitControl::metricIsPdf)
             .def("getMetric", &FitControl::getMetric);
             //.def_property("owner", &FitControl::getOwner, &FitControl::setOwner);

    py::class_<BinnedErrorFit,FitControl>(m, "BinnedErrorFit")
             .def(py::init<>());

}
