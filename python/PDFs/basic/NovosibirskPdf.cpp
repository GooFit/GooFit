#include <goofit/Python.h>

#include <goofit/PDFs/basic/NovosibirskPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/basic/NovosibirskPdf.h>

using namespace GooFit;

void init_NovosibirskPdf(py::module &m) {
    py::class_<NovosibirskPdf, GooPdf>(m, "NovosibirskPdf")
        .def(py::init<std::string, Observable, Variable, Variable, Variable>(),
             NovosibirskPdf_docs.c_str(),
             "name"_a,
             "x"_a,
             "m"_a,
             "s"_a,
             "t"_a)
        .def_static("help", []() { return HelpPrinter(NovosibirskPdf_docs); });
}
