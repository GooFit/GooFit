#include <goofit/Python.h>

#include <goofit/PDFs/basic/StepPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/basic/StepPdf.h>

using namespace GooFit;

void init_StepPdf(py::module &m) {
    py::class_<StepPdf, GooPdf>(m, "StepPdf")
        .def(py::init<std::string, Observable, Variable, int>(), StepPdf_docs.c_str(), "name"_a, "x"_a, "x0"_a, "Up"_a)
        .def_static("help", []() { return HelpPrinter(StepPdf_docs); });
}
