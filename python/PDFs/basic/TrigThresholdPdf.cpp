#include <goofit/Python.h>

#include <goofit/PDFs/basic/TrigThresholdPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/basic/TrigThresholdPdf.h>

using namespace GooFit;

void init_TrigThresholdPdf(py::module &m) {
    py::class_<TrigThresholdPdf, GooPdf>(m, "TrigThresholdPdf")
        .def(py::init<std::string, Observable, Variable, Variable, Variable, bool>(),
             TrigThresholdPdf_docs.c_str(),
             "name"_a,
             "x"_a,
             "thresh"_a,
             "trigConst"_a,
             "linConst"_a,
             "upper"_a = true)
        .def(py::init<std::string, Observable, Observable, Variable, Variable, Variable, Variable, bool>(),
             TrigThresholdPdf_docs.c_str(),
             "name"_a,
             "x"_a,
             "y"_a,
             "thresh"_a,
             "trigConst"_a,
             "linConst"_a,
             "massConstant"_a,
             "upper"_a)
        .def_static("help", []() { return HelpPrinter(TrigThresholdPdf_docs); });
}
