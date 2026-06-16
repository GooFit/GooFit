#include <pybind11/pybind11.h>

#include <pybind11/stl.h>

#include <Minuit2/FCNBase.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnUserParameters.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ROOT::Minuit2;

void init_MnMigrad(py::module &m) {
    // ROOT 6.34's Minuit2 dropped the (FCNBase, MnUserParameters, unsigned int) overload in favor of
    // (FCNBase, MnUserParameterState, MnStrategy); construct it explicitly so both APIs are supported.
    py::class_<MnMigrad, MnApplication>(m, "MnMigrad")
        .def(py::init([](const FCNBase &fcn, const MnUserParameters &par, unsigned int stra) {
                 return new MnMigrad(fcn, par, MnStrategy(stra));
             }),
             "fcn"_a,
             "par"_a,
             "stra"_a = 1);
}
