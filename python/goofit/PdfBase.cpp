#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/BinnedDataSet.h>
#include <goofit/PdfBase.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>

#include <iostream>

#include <Minuit2/FunctionMinimum.h>

namespace py = pybind11;
using namespace GooFit;
using namespace pybind11::literals;

void init_PdfBase(py::module &m) {
    py::class_<PdfBase>(m, "PdfBase")
        .def("setData", (void (PdfBase::*)(DataSet *)) & PdfBase::setData)
        .def("getParameters", &PdfBase::getParameters)
        .def("getParameterByName", &PdfBase::getParameterByName)
        .def("getObservables", &PdfBase::getObservables)
        .def("fitTo",
             &PdfBase::fitTo,
             py::call_guard<py::scoped_ostream_redirect>(),
             "Quick way to fit a PDF. Use a FitManager for more control.",
             "data"_a,
             "verbosity"_a = 3);
}
