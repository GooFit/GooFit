#include <pybind11/pybind11.h>

#include <goofit/Variable.h>
#include <goofit/BinnedDataSet.h>
#include <goofit/PDFs/basic/SmoothHistogramPdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_SmoothHistogramPdf(py::module &m) {
    py::class_<SmoothHistogramPdf, GooPdf>(m, "SmoothHistogramPdf")
        .def(py::init<std::string, BinnedDataSet *, Variable *>())
        ;
}








