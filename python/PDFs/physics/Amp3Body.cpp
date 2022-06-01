#include <goofit/Python.h>

#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <goofit/PDFs/physics/Amp3Body.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/physics/Amp3Body.h>

using namespace GooFit;

void init_Amp3Body(py::module &m) {
    py::class_<Amp3Body, Amp3BodyBase> cls(m, "Amp3Body");
    cls.def(py::init<std::string, Observable, Observable, EventNumber, DecayInfo3, GooPdf *>(),
            Amp3Body_docs.c_str(),
            "name"_a,
            "m12"_a,
            "m13"_a,
            "eventNumber"_a,
            "decay"_a,
            "eff"_a,
            py::keep_alive<1, 6>(), // Important to keep decay alive, to keep PDFs alive
            py::keep_alive<1, 7>())
        .def("setDataSize", &Amp3Body::setDataSize, "dataSize"_a, "evtSize"_a = 3, "offset"_a = 0)
        .def("getCachedWave", &Amp3Body::getCachedWave, "i"_a)
        .def(
            "sumCachedWave",
            [](Amp3Body &self, size_t i) { return std::complex<fptype>(self.sumCachedWave(i)); },
            "i"_a)
        .def("fit_fractions",
             &Amp3Body::fit_fractions,
             "Using the current dataset, return the cached fit fraction values")
        .def("getDecayInfo", &Amp3Body::getDecayInfo, "Return DecayInfo")
        .def_static("resetCacheCounter", &Amp3Body::resetCacheCounter)
        .def("normalize", &Amp3Body::normalize)
        .def_static("help", []() { return HelpPrinter(Amp3Body_docs); })
        .def("setGenerationOffset", &Amp3Body::setGenerationOffset, "off"_a)
        .def("getGenerationOffset", &Amp3Body::getGenerationOffset, pybind11::return_value_policy::copy)
        .def("GenerateSig", [](Amp3Body &self, size_t numEvents) {
            mcbooster::ParticlesSet_h particles; // typedef for std::vector<Particles_h *>,
            mcbooster::VariableSet_h variables;
            mcbooster::RealVector_h weights;
            mcbooster::BoolVector_h flags;

            std::tie(particles, variables, weights, flags) = self.GenerateSig(numEvents);

            py::array_t<fptype> pyparticles{{(size_t)3, numEvents, (size_t)3}};
            py::array_t<fptype> pyvariables{{(size_t)3, numEvents}};
            py::array_t<fptype> pyweights{static_cast<ssize_t>(numEvents)};
            py::array_t<bool> pyflags{static_cast<ssize_t>(numEvents)};

            for(int i = 0; i < 3; i++) {
                for(int j = 0; j < numEvents; j++) {
                    for(int k = 0; k < 3; k++) {
                        pyparticles.mutable_at(i, j, k) = (*(particles[i]))[j].get(k);
                    }
                }
            }

            for(int i = 0; i < 3; i++) {
                for(int j = 0; j < numEvents; j++) {
                    pyvariables.mutable_at(i, j) = (*(variables[i]))[j];
                }
            }

            for(int i = 0; i < numEvents; i++) {
                pyweights.mutable_at(i) = weights[i];
            }

            for(int i = 0; i < numEvents; i++) {
                pyflags.mutable_at(i) = flags[i];
            }

            return std::make_tuple(pyparticles, pyvariables, pyweights, pyflags);
        });

    m.attr("DalitzPlotPdf") = cls;
}
