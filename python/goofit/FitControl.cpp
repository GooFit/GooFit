#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <iostream>
#include <memory>

#include <goofit/FitControl.h>

#include <Minuit2/FunctionMinimum.h>

namespace py = pybind11;
using namespace GooFit;
using namespace pybind11::literals;

void init_FitControl(py::module &m) {
    py::class_<FitControl, std::shared_ptr<FitControl>>(m, "FitControl")
        .def("binnedFit", &FitControl::binnedFit)
        .def("binErrors", &FitControl::binErrors)
        .def("metricIsPdf", &FitControl::metricIsPdf)
        .def("getMetric", &FitControl::getMetric);

    py::class_<BinnedErrorFit, FitControl, std::shared_ptr<BinnedErrorFit>>(m, "BinnedErrorFit").def(py::init<>());
    py::class_<UnbinnedNllFit, FitControl, std::shared_ptr<UnbinnedNllFit>>(m, "UnbinnedNllFit").def(py::init<>());
    py::class_<BinnedNllFit, FitControl, std::shared_ptr<BinnedNllFit>>(m, "BinnedNllFit").def(py::init<>());
    py::class_<BinnedChisqFit, FitControl, std::shared_ptr<BinnedChisqFit>>(m, "BinnedChisqFit").def(py::init<>());
}
