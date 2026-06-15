#include <pybind11/pybind11.h>

#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <Minuit2/MnUserCovariance.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ROOT::Minuit2;

void init_MnUserCovariance(py::module &m) {
    py::class_<MnUserCovariance>(m, "MnUserCovariance")
        .def("Nrow", &MnUserCovariance::Nrow)
        .def("size", &MnUserCovariance::size)
        .def("__call__",
             (double (MnUserCovariance::*)(unsigned int, unsigned int) const) & MnUserCovariance::operator())
        .def("to_matrix", [](const MnUserCovariance &self) {
            size_t n = self.Nrow();
            py::array_t<double> matrix({n, n});
            for(size_t i = 0; i < n; i++)
                for(size_t j = 0; j < n; j++)
                    matrix.mutable_at(i, j) = self(i, j);
            return matrix;
        });
}
