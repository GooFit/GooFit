#include <goofit/Python.h>

#include <goofit/PDFs/basic/LandauPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/basic/LandauPdf.h>

using namespace GooFit;

void init_LandauPdf(py::module &m) {
    py::class_<LandauPdf, GooPdf>(m, "LandauPdf")
        .def(py::init<std::string, Observable, Variable, Variable>(),
             LandauPdf_docs.c_str(),
             "name"_a,
             "x"_a,
             "mpv"_a,
             "sigma"_a)
        .def_static("help", []() { return HelpPrinter(LandauPdf_docs); });
}
