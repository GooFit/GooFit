#include <goofit/Python.h>

#include <pybind11/stl.h>

#include <goofit/PDFs/basic/InterHistPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/basic/InterHistPdf.h>

using namespace GooFit;

void init_InterHistPdf(py::module &m) {
    py::class_<InterHistPdf, GooPdf>(m, "InterHistPdf")
        .def(py::init<std::string, BinnedDataSet *, std::vector<Variable>, std::vector<Observable>>(),
             InterHistPdf_docs.c_str(),
             "name"_a,
             "x"_a,
             "params"_a,
             "obses"_a,
             py::keep_alive<1, 3>())
        .def_static("help", []() { return HelpPrinter(InterHistPdf_docs); });
}
