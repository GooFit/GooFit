#include <goofit/Python.h>

#include <goofit/PDFs/basic/BWPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/basic/BWPdf.h>

using namespace GooFit;

void init_BWPdf(py::module &m) {
    py::class_<BWPdf, GooPdf>(m, "BWPdf")
        .def(py::init<std::string, Observable, Variable, Variable>(),
             "A non-relativistic Breit-Wigner function, sometimes called "
             "a Cauchy function",
             "name"_a,
             "x"_a,
             "mean"_a,
             "Gamma"_a)
        .def_static("help", []() { return HelpPrinter(BWPdf_docs); });
}
