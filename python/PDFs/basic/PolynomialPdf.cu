#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/PDFs/basic/PolynomialPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;
using namespace pybind11::literals;

void init_PolynomialPdf(py::module &m) {
    py::class_<PolynomialPdf, GooPdf>(m, "PolynomialPdf")
        .def(py::init<std::string, Variable *, std::vector<Variable *>, Variable *, unsigned int>(),
             "n"_a,
             "_x"_a,
             "weights"_a,
             "x0"_a           = nullptr,
             "lowestDegree"_a = 0)
        .def(py::init<std::string,
                      std::vector<Variable *>,
                      std::vector<Variable *>,
                      std::vector<Variable *>,
                      unsigned int>(),
             "n"_a,
             "obses"_a,
             "coeffs"_a,
             "offsets"_a,
             "maxDegree"_a);
}
