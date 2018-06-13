#include <goofit/Python.h>

#include <pybind11/stl.h>

#include <goofit/PDFs/combine/ProdPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/combine/ProdPdf.h>

using namespace GooFit;

void init_ProdPdf(py::module &m) {
    py::class_<ProdPdf, CombinePdf>(m, "ProdPdf")
        .def(py::init<std::string, std::vector<PdfBase *>>(), ProdPdf_docs.c_str(), "name"_a, "comps"_a)

        .def_static("help", []() { return HelpPrinter(ProdPdf_docs); })

        ;
}
