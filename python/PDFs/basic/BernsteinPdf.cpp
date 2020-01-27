#include <goofit/Python.h>

#include <pybind11/stl.h>

#include <goofit/PDFs/basic/BernsteinPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/basic/BernsteinPdf.h>

using namespace GooFit;

void init_BernsteinPdf(py::module &m) {
    py::class_<BernsteinPdf, GooPdf>(m, "BernsteinPdf")
        .def(py::init<std::string, Observable, std::vector<Variable>, int>(),
             "Bernstein form, that is a linear combination of "
             "Bernstein basis polynomials",
             "name"_a,
             "obs"_a,
             "coeffs"_a,
             "mindeg"_a)
        .def_static("help", []() { return HelpPrinter(BernsteinPdf_docs); });
}
