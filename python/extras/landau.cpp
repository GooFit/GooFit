#include <pybind11/pybind11.h>

#include <pybind11/numpy.h>

#include <goofit/cpp/landau.h>

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(landau, m) {
    m.def("landauPDF", py::vectorize(landauPDF), "x"_a, "x0"_a = 0., "xi"_a = 1.);
    m.def("gaussPDF", py::vectorize(gaussPDF), "x"_a, "mu"_a = 0., "sigma"_a = 1.);
    m.def("landauGaussPDF", py::vectorize(landauGaussPDF), "x"_a, "mu"_a = 0., "eta"_a = 1., "sigma"_a = 1.);
    m.def("landau_quantile", py::vectorize(landau_quantile), "z"_a, "xi"_a = 1.0);
    m.def("landau_quantile_c", py::vectorize(landau_quantile_c), "z"_a, "xi"_a = 1.0);
}
