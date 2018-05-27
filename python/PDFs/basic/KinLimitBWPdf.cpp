#include <goofit/Python.h>

#include <goofit/PDFs/basic/KinLimitBWPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/basic/KinLimitBWPdf.h>

using namespace GooFit;

void init_KinLimitBWPdf(py::module &m) {
    py::class_<KinLimitBWPdf, GooPdf>(m, "KinLimitBWPdf")
        .def(py::init<std::string, Observable, Variable, Variable>(),
             KinLimitBWPdf_docs.c_str(),
             "name"_a,
             "x"_a,
             "m"_a,
             "s"_a)
        .def_static("help", []() { return HelpPrinter(KinLimitBWPdf_docs); });
}
