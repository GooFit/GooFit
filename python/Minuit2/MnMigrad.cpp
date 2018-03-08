#include <pybind11/pybind11.h>

#include <pybind11/stl.h>

#include <Minuit2/FCNBase.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnUserParameters.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ROOT::Minuit2;

void init_MnMigrad(py::module &m) {
    py::class_<MnMigrad, MnApplication>(m, "MnMigrad")
        .def(py::init<const FCNBase &, const MnUserParameters &, unsigned int>(), "fcn"_a, "par"_a, "stra"_a = 1);
}
