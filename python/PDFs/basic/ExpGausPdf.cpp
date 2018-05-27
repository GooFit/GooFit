#include <goofit/Python.h>

#include <goofit/PDFs/basic/ExpGausPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/basic/ExpGausPdf.h>

using namespace GooFit;

void init_ExpGausPdf(py::module &m) {
    py::class_<ExpGausPdf, GooPdf>(m, "ExpGausPdf")
        .def(py::init<std::string, Observable, Variable, Variable, Variable>(),
             "An exponential decay convolved with a Gaussian resolution.",
             "name"_a,
             "x"_a,
             "m"_a,
             "s"_a,
             "t"_a)
        .def_static("help", []() { return HelpPrinter(ExpGausPdf_docs); });
    ;
}
