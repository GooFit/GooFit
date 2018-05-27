#include <goofit/Python.h>

#include <goofit/BinnedDataSet.h>
#include <goofit/PDFs/basic/SmoothHistogramPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/basic/SmoothHistogramPdf.h>

using namespace GooFit;

void init_SmoothHistogramPdf(py::module &m) {
    py::class_<SmoothHistogramPdf, GooPdf>(m, "SmoothHistogramPdf")
        .def(py::init<std::string, BinnedDataSet *, Variable>(),
             SmoothHistogramPdf_docs.c_str(),
             "name"_a,
             "x"_a,
             "smoothing"_a,
             py::keep_alive<1, 3>())
        .def_static("help", []() { return HelpPrinter(SmoothHistogramPdf_docs); });
}
