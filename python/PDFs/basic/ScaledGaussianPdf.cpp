#include <goofit/Python.h>

#include <goofit/PDFs/basic/ScaledGaussianPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/basic/ScaledGaussianPdf.h>

using namespace GooFit;

void init_ScaledGaussianPdf(py::module &m) {
    py::class_<ScaledGaussianPdf, GooPdf>(m, "ScaledGaussianPdf")
        .def(py::init<std::string, Observable, Variable, Variable, Variable, Variable>(),
             ScaledGaussianPdf_docs.c_str(),
             "name"_a,
             "x"_a,
             "m"_a,
             "s"_a,
             "d"_a,
             "e"_a)
        .def_static("help", []() { return HelpPrinter(ScaledGaussianPdf_docs); });
}
