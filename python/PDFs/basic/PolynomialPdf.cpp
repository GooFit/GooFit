#include <goofit/Python.h>

#include <pybind11/stl.h>

#include <goofit/PDFs/basic/PolynomialPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/basic/PolynomialPdf.h>

using namespace GooFit;

void init_PolynomialPdf(py::module &m) {
    py::class_<PolynomialPdf, GooPdf>(m, "PolynomialPdf")
        .def(py::init<std::string, Observable, std::vector<Variable>, Variable, unsigned int>(),
             PolynomialPdf_docs.c_str(),
             "name"_a,
             "x"_a,
             "weights"_a,
             "x0"_a,
             "lowestDegree"_a = 0)
        .def(py::init<std::string, Observable, std::vector<Variable>, unsigned int>(),
             PolynomialPdf_docs.c_str(),
             "name"_a,
             "x"_a,
             "weights"_a,
             "lowestDegree"_a = 0)
        .def(py::init<std::string,
                      std::vector<Observable>,
                      std::vector<Variable>,
                      std::vector<Variable>,
                      unsigned int>(),
             PolynomialPdf_docs.c_str(),
             "name"_a,
             "obses"_a,
             "coeffs"_a,
             "offsets"_a,
             "maxDegree"_a = 0)
        .def_static("help", []() { return HelpPrinter(PolynomialPdf_docs); });
}
