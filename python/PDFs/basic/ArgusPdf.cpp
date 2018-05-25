#include <goofit/Python.h>

#include <goofit/PDFs/basic/ArgusPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/basic/ArgusPdf.h>

using namespace GooFit;

void init_ArgusPdf(py::module &m) {
    py::class_<ArgusPdf, GooPdf>(m, "ArgusPdf")
        .def(py::init<std::string, Observable, Variable, Variable, bool>(),
             "Standard Argus",
             "n"_a,
             "x"_a,
             "m"_a,
             "s"_a,
             "upper"_a)

        .def(py::init<std::string, Observable, Variable, Variable, bool, Variable>(),
             "Power version of Argus",
             "n"_a,
             "x"_a,
             "m"_a,
             "s"_a,
             "upper"_a,
             "power"_a)

        .def_static("help", []() { return HelpPrinter(ArgusPdf_docs); });
}
