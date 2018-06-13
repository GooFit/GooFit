#include <goofit/Python.h>

#include <goofit/PDFs/combine/CompositePdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/combine/CompositePdf.h>

using namespace GooFit;

void init_CompositePdf(py::module &m) {
    py::class_<CompositePdf, CombinePdf>(m, "CompositePdf")
        .def(py::init<std::string, PdfBase *, PdfBase *>(),
             CompositePdf_docs.c_str(),
             "name"_a,
             "core"_a,
             "shell"_a,
             py::keep_alive<1, 3>(),
             py::keep_alive<1, 4>())

        .def_static("help", []() { return HelpPrinter(CompositePdf_docs); })

        ;
}
