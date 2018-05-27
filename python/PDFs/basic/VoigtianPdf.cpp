#include <goofit/Python.h>

#include <goofit/PDFs/basic/VoigtianPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/basic/VoigtianPdf.h>

using namespace GooFit;

void init_VoigtianPdf(py::module &m) {
    py::class_<VoigtianPdf, GooPdf>(m, "VoigtianPdf")
        .def(py::init<std::string, Observable, Variable, Variable, Variable>(),
             VoigtianPdf_docs.c_str(),
             "name"_a,
             "x"_a,
             "m"_a,
             "s"_a,
             "w"_a)
        .def_static("help", []() { return HelpPrinter(VoigtianPdf_docs); });
}
