#include <goofit/Python.h>

#include <goofit/PDFs/basic/JohnsonSUPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/basic/JohnsonSUPdf.h>

using namespace GooFit;

void init_JohnsonSUPdf(py::module &m) {
    py::class_<JohnsonSUPdf, GooPdf>(m, "JohnsonSUPdf")
        .def(py::init<std::string, Observable, Variable, Variable, Variable, Variable>(),
             JohnsonSUPdf_docs.c_str(),
             "name"_a,
             "x"_a,
             "m"_a,
             "s"_a,
             "g"_a,
             "d"_a)
        .def_static("help", []() { return HelpPrinter(JohnsonSUPdf_docs); });
}
