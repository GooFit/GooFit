#include <goofit/Python.h>

#include <pybind11/stl.h>

#include <goofit/PDFs/combine/EventWeightedAddPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/combine/EventWeightedAddPdf.h>

using namespace GooFit;

void init_EventWeightedAddPdf(py::module &m) {
    py::class_<EventWeightedAddPdf, CombinePdf>(m, "EventWeightedAddPdf")
        .def(py::init<std::string, std::vector<Observable>, std::vector<PdfBase *>>(),
             EventWeightedAddPdf_docs.c_str(),
             "name"_a,
             "weights"_a,
             "comps"_a)

        .def_static("help", []() { return HelpPrinter(EventWeightedAddPdf_docs); })

        ;
}
