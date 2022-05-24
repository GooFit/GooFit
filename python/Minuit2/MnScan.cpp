#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <sstream>

#include <Minuit2/MnScan.h>
#include <Minuit2/MnPrint.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ROOT::Minuit2;

void init_MnScan(py::module &m) {
    py::class_<MnScan>(m, "MnScan")
        .def("Scan", &MnScan::Scan, "Scan parameter.", "par"_a, "maxsteps"_a = 41, "low"_a = 0, "high"_a = 0);
}
