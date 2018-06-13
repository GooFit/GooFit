#include <goofit/Python.h>

#include <pybind11/stl.h>

#include <goofit/PDFs/combine/AddPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/combine/AddPdf.h>

using namespace GooFit;

void init_AddPdf(py::module &m) {
    py::class_<AddPdf, CombinePdf>(m, "AddPdf")
        .def(py::init<std::string, std::vector<Variable>, std::vector<PdfBase *>>(),
             AddPdf_docs.c_str(),
             "name"_a,
             "weights"_a,
             "comps"_a)

        .def(py::init<std::string, Variable, PdfBase *, PdfBase *>(),
             AddPdf_docs.c_str(),
             "name"_a,
             "frac1"_a,
             "func1"_a,
             "func2"_a,
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>())

        .def_static("help", []() { return HelpPrinter(AddPdf_docs); })

        ;
}
