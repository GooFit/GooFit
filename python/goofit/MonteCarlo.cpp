#include <pybind11/pybind11.h>

#include <goofit/utilities/MonteCarlo.h>

namespace py = pybind11;
using namespace fmt::literals;
using namespace GooFit;

void init_MonteCarlo(py::module &m) {
    m.def("fillDataSetMC1D", &fillDataSetMC1D, "Fill a dataset with 1D Monte Carlo");
    //            "pdf"_a, "var"_a, "nTotal"_a, "seed"_a = 0);
}
