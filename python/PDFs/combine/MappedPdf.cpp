#include <goofit/Python.h>

#include <pybind11/stl.h>

#include <goofit/PDFs/combine/MappedPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/combine/MappedPdf.h>

using namespace GooFit;

void init_MappedPdf(py::module &m) {
    py::class_<MappedPdf, CombinePdf>(m, "MappedPdf")
        .def(py::init<std::string, GooPdf *, std::vector<GooPdf *> &>(),
             MappedPdf_docs.c_str(),
             "name"_a,
             "m"_a,
             "t"_a,
             py::keep_alive<1, 3>(),
             py::keep_alive<1, 4>())

        .def_static("help", []() { return HelpPrinter(MappedPdf_docs); })

        ;
}
