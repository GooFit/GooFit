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

    py::class_<DecayInfo3>(m, "DecayInfo3")
        .def(py::init<>())
        .def_readwrite("motherMass", &DecayInfo3::motherMass)
        .def_readwrite("daug1Mass", &DecayInfo3::daug1Mass)
        .def_readwrite("daug2Mass", &DecayInfo3::daug2Mass)
        .def_readwrite("daug3Mass", &DecayInfo3::daug3Mass)
        .def_readwrite("meson_radius", &DecayInfo3::meson_radius)
        .def_property("resonances",
                      [](DecayInfo3 &self) { return self.resonances; },
                      py::cpp_function([](DecayInfo3 &self, std::vector<ResonancePdf *> val) { self.resonances = val; },
                                       py::keep_alive<1, 2>()))
        .def("add_resonance",
             [](DecayInfo3 &self, ResonancePdf *toadd) { self.resonances.push_back(toadd); },
             "Append a resonance",
             "resonance"_a,
             py::keep_alive<1, 2>())
        .def("__str__", [](DecayInfo3 &self) {
            return fmt::format("M={} GeV\nm1={}\nGeV\nm2={}\nm3={}\nrad={}\nN res: {}",
                               self.motherMass,
                               self.daug1Mass,
                               self.daug2Mass,
                               self.daug3Mass,
                               self.meson_radius,
                               self.resonances.size());
        });

    py::class_<DecayInfo3t, DecayInfo3>(m, "DecayInfo3t")
        .def(py::init<Variable, Variable, Variable>(), "tau"_a, "xmixing"_a, "ymixing"_a)
        .def_readonly("_tau", &DecayInfo3t::_tau)
        .def_readonly("_xmixing", &DecayInfo3t::_xmixing)
        .def_readonly("_ymixing", &DecayInfo3t::_ymixing);

    py::class_<DecayInfo4>(m, "DecayInfo4")
        .def(py::init<>())
        .def_readwrite("particle_masses", &DecayInfo4::particle_masses)
        .def_readwrite("meson_radius", &DecayInfo4::meson_radius)
        .def_readwrite("amplitudes", &DecayInfo4::amplitudes)
        .def_readwrite("amplitudes_B", &DecayInfo4::amplitudes_B)
        .def_property("amplitudes",
                      [](DecayInfo4 &self) { return self.amplitudes; },
                      py::cpp_function([](DecayInfo4 &self, std::vector<Amplitude *> val) { self.amplitudes = val; },
                                       py::keep_alive<1, 2>()))
        .def_property("amplitudes_B",
                      [](DecayInfo4 &self) { return self.amplitudes_B; },
                      py::cpp_function([](DecayInfo4 &self, std::vector<Amplitude *> val) { self.amplitudes_B = val; },
                                       py::keep_alive<1, 2>()));

    py::class_<DecayInfo4t, DecayInfo4>(m, "DecayInfo4t")
        .def(py::init<Variable, Variable, Variable, Variable>(), "tau"_a, "xmixing"_a, "ymixing"_a, "SqWStoRSrate"_a)
        .def_readonly("_tau", &DecayInfo4t::_tau)
        .def_readonly("_xmixing", &DecayInfo4t::_xmixing)
        .def_readonly("_ymixing", &DecayInfo4t::_ymixing)
        .def_readonly("_SqWStoRSrate", &DecayInfo4t::_SqWStoRSrate);
}
