#include <goofit/Python.h>

#include <goofit/PDFs/basic/BifurGaussPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/basic/BifurGaussPdf.h>

using namespace GooFit;

void init_BifurGaussPdf(py::module &m) {
    py::class_<BifurGaussPdf, GooPdf>(m, "BifurGaussPdf")
        .def(py::init<std::string, Observable, Variable, Variable, Variable>(),
             "A two-sided Gaussian, with a sigma that varies "
             "depending on which side of the mean you are on",
             "name"_a,
             "x"_a,
             "m"_a,
             "sL"_a,
             "sR"_a)
        .def_static("help", []() { return HelpPrinter(BifurGaussPdf_docs); });
}
