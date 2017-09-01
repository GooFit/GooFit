#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;
using namespace pybind11::literals;

void init_DalitzPlotHelpers(py::module &m) {
    py::enum_<DaughterPair>(m, "DaughterPair")
        .value("PAIR_12", DaughterPair::PAIR_12)
        .value("PAIR_12", DaughterPair::PAIR_13)
        .value("PAIR_23", DaughterPair::PAIR_23)
        .export_values();

    py::class_<DecayInfo>(m, "DecayInfo")
        .def_readwrite("motherMass", &DecayInfo::motherMass)
        .def_readwrite("daug1Mass", &DecayInfo::daug1Mass)
        .def_readwrite("daug2Mass", &DecayInfo::daug2Mass)
        .def_readwrite("daug3Mass", &DecayInfo::daug3Mass)
        .def_readwrite("meson_radius", &DecayInfo::meson_radius)
        .def_readwrite("_tau", &DecayInfo::_tau)
        .def_readwrite("_xmixing", &DecayInfo::_xmixing)
        .def_readwrite("_ymixing", &DecayInfo::_ymixing)
        .def("add_resonance",
             [](DecayInfo &self, ResonancePdf *toadd) { self.resonances.push_back(toadd); },
             "Append a resonance",
             "resonance"_a);
}
