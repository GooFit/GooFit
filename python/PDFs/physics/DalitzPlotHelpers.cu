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
            .value("PAIR_13", DaughterPair::PAIR_13)
            .value("PAIR_23", DaughterPair::PAIR_23)
            .export_values();

    py::enum_<DP4Pair>(m, "DP4Pair")
            .value("M_12", DP4Pair::M_12)
            .value("M_34", DP4Pair::M_34)
            .value("M_13", DP4Pair::M_13)
            .value("M_14", DP4Pair::M_14)
            .value("M_23", DP4Pair::M_23)
            .value("M_24", DP4Pair::M_24)
            .value("M_12_3", DP4Pair::M_12_3)
            .value("M_13_2", DP4Pair::M_13_2)
            .value("M_23_1", DP4Pair::M_23_1)
            .value("M_12_4", DP4Pair::M_12_4)
            .value("M_14_2", DP4Pair::M_14_2)
            .value("M_24_1", DP4Pair::M_24_1)
            .value("M_13_4", DP4Pair::M_13_4)
            .value("M_14_3", DP4Pair::M_14_3)
            .value("M_34_1", DP4Pair::M_34_1)
            .value("M_23_4", DP4Pair::M_23_4)
            .value("M_24_3", DP4Pair::M_24_3)
            .value("M_34_2", DP4Pair::M_34_2)
            .export_values();

    py::class_<DecayInfo>(m, "DecayInfo")
            .def(py::init<>())
            .def_readwrite("motherMass", &DecayInfo::motherMass)
            .def_readwrite("daug1Mass", &DecayInfo::daug1Mass)
            .def_readwrite("daug2Mass", &DecayInfo::daug2Mass)
            .def_readwrite("daug3Mass", &DecayInfo::daug3Mass)
            .def_readwrite("meson_radius", &DecayInfo::meson_radius)
            .def_readwrite("_tau", &DecayInfo::_tau)
            .def_readwrite("_xmixing", &DecayInfo::_xmixing)
            .def_readwrite("_ymixing", &DecayInfo::_ymixing)
            .def_property("resonances", [](DecayInfo &self){return self.resonances;},
                py::cpp_function([](DecayInfo &self, std::vector<ResonancePdf*> val){self.resonances = val;},
                py::keep_alive<1,2>()))
            .def("add_resonance",
                 [](DecayInfo &self, ResonancePdf *toadd) { self.resonances.push_back(toadd); },
                        "Append a resonance",
                        "resonance"_a);

    py::class_<DecayInfo_DP>(m, "DecayInfo_DP")
            .def(py::init<>())
            .def_readwrite("particle_masses", &DecayInfo_DP::particle_masses)
            .def_readwrite("meson_radius", &DecayInfo_DP::meson_radius)
            .def_readwrite("amplitudes", &DecayInfo_DP::amplitudes)
            .def_readwrite("amplitudes_B", &DecayInfo_DP::amplitudes_B)
            .def_readwrite("_tau", &DecayInfo_DP::_tau)
            .def_readwrite("_xmixing", &DecayInfo_DP::_xmixing)
            .def_readwrite("_ymixing", &DecayInfo_DP::_ymixing)
            .def_readwrite("_SqWStoRSrate", &DecayInfo_DP::_SqWStoRSrate)
            .def_property("amplitudes", [](DecayInfo_DP &self){return self.amplitudes;},
                py::cpp_function([](DecayInfo_DP &self, std::vector<Amplitude*> val){self.amplitudes = val;},
                py::keep_alive<1,2>()))
            .def_property("amplitudes_B", [](DecayInfo_DP &self){return self.amplitudes_B  ;},
                py::cpp_function([](DecayInfo_DP &self, std::vector<Amplitude*> val){self.amplitudes_B = val;},
                py::keep_alive<1,2>()));
}
