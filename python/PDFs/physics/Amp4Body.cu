#include <goofit/Python.h>

#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <goofit/DataSet.h>
#include <goofit/PDFs/physics/DP4Pdf.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/Variable.h>

using namespace GooFit;

void init_Amp4Body(py::module &m) {
    py::class_<Amp4Body, Amp4BodyBase> cls(m, "Amp4Body");
    cls.def(py::init<std::string, std::vector<Observable>, DecayInfo4, GooPdf *>(),
            "name"_a,
            "observables"_a,
            "decay"_a,
            "eff"_a,
            py::keep_alive<1, 4>(),
            py::keep_alive<1, 5>())
        .def(py::init<std::string, std::vector<Observable>, DecayInfo4, GooPdf *, unsigned int>(),
             "name"_a,
             "observables"_a,
             "decay"_a,
             "eff"_a,
             "MCeventsNorm"_a,
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>())
        .def("normalize", &DPPdf::normalize)
        .def("setDataSize", &DPPdf::setDataSize, "dataSize"_a, "evtSize"_a = 6)
        .def("setForceIntegrals", &DPPdf::setForceIntegrals, "fset"_a = true)
        .def("getMCevents", &DPPdf::getMCevents)
        .def("setGenerationOffset", &DPPdf::setGenerationOffset, "off"_a)
        .def("GenerateSig",
             [](DPPdf &self, size_t numEvents) {
                 mcbooster::ParticlesSet_h particles; // typedef for std::vector<Particles_h *>,
                 mcbooster::VariableSet_h variables;
                 mcbooster::RealVector_h weights;
                 mcbooster::BoolVector_h flags;

                 std::tie(particles, variables, weights, flags) = self.GenerateSig(numEvents);

                 py::array_t<fptype> pyparticles{{(size_t)4, numEvents, (size_t)4}};
                 py::array_t<fptype> pyvariables{{(size_t)5, numEvents}};
                 py::array_t<fptype> pyweights{(py::ssize_t)numEvents};
                 py::array_t<bool> pyflags{(py::ssize_t)numEvents};

                 for(int i = 0; i < 4; i++) {
                     for(int j = 0; j < numEvents; j++) {
                         for(int k = 0; k < 4; k++) {
                             pyparticles.mutable_at(i, j, k) = (*(particles[i]))[j].get(k);
                         }
                     }
                 }

                 for(int i = 0; i < 5; i++) {
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
             })

        ;

    m.attr("DPPdf") = cls;
}
