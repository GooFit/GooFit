#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/basic/PolynomialPdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_PolynomialPdf(py::module &m) {
    py::class_<PolynomialPdf, GooPdf>(m, "PolynomialPdf")
        .def(py::init<std::string, Variable *, std::vector<Variable *>>())
        .def(py::init<std::string, Variable *, std::vector<Variable *>, Variable *>())

        .def(py::init<std::string, Variable *, std::vector<Variable *>, Variable *, unsigned int >())
        .def(py::init<std::string, std::vector<Variable *>, std::vector<Variable *>, std::vector<Variable *>, unsigned int>())
        ;
}
