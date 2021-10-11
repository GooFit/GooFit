#include <goofit/Python.h>

#include <pybind11/iostream.h>
#include <pybind11/stl.h>

#include <goofit/BinnedDataSet.h>
#include <goofit/PdfBase.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>

#include <iostream>

#include <Minuit2/FunctionMinimum.h>

using namespace GooFit;

void init_PdfBase(py::module &m) {
    py::class_<PdfBase>(m, "PdfBase")
        .def("setData", (void(PdfBase::*)(DataSet *)) & PdfBase::setData, "ptr"_a = nullptr)
        .def("getData", &PdfBase::getData)
        .def("getName", &PdfBase::getName)
        .def("getParameters", &PdfBase::getParameters)
        .def("getParameterByName", &PdfBase::getParameterByName)
        .def("getObservables", &PdfBase::getObservables)
        .def("getFunctionIndex", &PdfBase::getFunctionIndex)
        .def("__str__",
             [](const PdfBase &pdf) {
                 std::stringstream str;
                 str << pdf;
                 return str.str();
             })
        .def("fillMCDataSimple",
             &PdfBase::fillMCDataSimple,
             "Fill in events in the dataset. Must have a dataset and"
             " it will be appended to.",
             "events"_a,
             "seed"_a = 0)
        .def("fitTo",
             &PdfBase::fitTo,
             py::call_guard<py::scoped_ostream_redirect>(),
             "Quick way to fit a PDF. Use a FitManager for more control.",
             "data"_a,
             "verbosity"_a = 3)
        .def("fit",
             &PdfBase::fit,
             py::call_guard<py::scoped_ostream_redirect>(),
             "Fit to an existing attached dataset.",
             "verbosity"_a = 3);
}
