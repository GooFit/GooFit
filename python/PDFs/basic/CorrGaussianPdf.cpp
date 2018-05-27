#include <goofit/Python.h>

#include <goofit/PDFs/basic/CorrGaussianPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/basic/CorrGaussianPdf.h>

using namespace GooFit;

void init_CorrGaussianPdf(py::module &m) {
    py::class_<CorrGaussianPdf, GooPdf>(m, "CorrGaussianPdf")
        .def(py::init<std::string, Observable, Observable, Variable, Variable, Variable, Variable, Variable>(),
             "A correlated Gaussian - that is, a function of "
             "two variables x and y, each described by a Gaussian "
             "distribution, but the width of the y distribution depends on x.",
             "name"_a,
             "x"_a,
             "y"_a,
             "mean1"_a,
             "sigma1"_a,
             "mean2"_a,
             "sigma2"_a,
             "correlation"_a)
        .def_static("help", []() { return HelpPrinter(CorrGaussianPdf_docs); });
}
