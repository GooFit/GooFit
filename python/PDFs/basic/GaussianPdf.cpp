#include <goofit/Python.h>

#include <goofit/PDFs/basic/GaussianPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/basic/GaussianPdf.h>

using namespace GooFit;

void init_GaussianPdf(py::module &m) {
    py::class_<GaussianPdf, GooPdf>(m, "GaussianPdf")
        .def(py::init<std::string, Observable, Variable, Variable>(),
             GaussianPdf_docs.c_str(),
             "name"_a,
             "x"_a,
             "m"_a,
             "s"_a)
        .def_static("help", []() { return HelpPrinter(GaussianPdf_docs); });
}
