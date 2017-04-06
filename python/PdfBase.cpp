
#include <pybind11/pybind11.h>
#include <goofit/PdfBase.h>

namespace py = pybind11;

void init_PdfBase(py::module &m) {
    py::class_<PdfBase>(m, "PdfBase")
        .def("setData", (void (PdfBase::*)(UnbinnedDataSet*)) &PdfBase::setData)
        .def("setData", (void (PdfBase::*)(BinnedDataSet*)) &PdfBase::setData)
        .def("setData", (void (PdfBase::*)(std::vector<std::map<Variable*, fptype>>&)) &PdfBase::setData)
        ;

}

