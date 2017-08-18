
#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>
#include <iostream>

#include <goofit/Variable.h>
#include <goofit/BinnedDataSet.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/PdfBase.h>

#include <Minuit2/FunctionMinimum.h>

namespace py = pybind11;
using namespace GooFit;
using namespace pybind11::literals;

void init_PdfBase(py::module &m) {
    py::class_<PdfBase>(m, "PdfBase")
        .def("setData", (void (PdfBase::*)(DataSet *)) & PdfBase::setData)
        .def("setData", (void (PdfBase::*)(std::vector<std::map<Variable *, fptype>> &)) & PdfBase::setData)
        //.def("fitTo", &PdfBase::fitTo) <- add Minuit bindings to make this work
        .def("fitTo", [](PdfBase& self, DataSet* data, int verbosity){
                py::scoped_output_redirect redir(
                    std::cout, py::module::import("sys").attr("stdout")
                );
                self.fitTo(data, verbosity);
                return;
            },
            "Quick way to fit a PDF. Use a FitManager for more control.", 
            "data"_a, "verbosity"_a = 3)
    ;
}
