#include <goofit/Python.h>

#include <goofit/PDFs/combine/ConvolutionPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/combine/ConvolutionPdf.h>

using namespace GooFit;

void init_ConvolutionPdf(py::module &m) {
    py::class_<ConvolutionPdf, CombinePdf>(m, "ConvolutionPdf")
        .def(py::init<std::string, Observable, GooPdf *, GooPdf *>(),
             ConvolutionPdf_docs.c_str(),
             "name"_a,
             "x"_a,
             "model"_a,
             "resolution"_a,
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>())
        .def(py::init<std::string, Observable, GooPdf *, GooPdf *, unsigned int>(),
             ConvolutionPdf_docs.c_str(),
             "name"_a,
             "x"_a,
             "model"_a,
             "resolution"_a,
             "numOthers"_a,
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>())

        .def("setIntegrationConstants", &ConvolutionPdf::setIntegrationConstants, "lo", "hi", "step")
        .def_static("help", []() { return HelpPrinter(ConvolutionPdf_docs); })

        ;
}
