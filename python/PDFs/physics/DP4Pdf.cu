#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/DataSet.h>
#include <goofit/PDFs/physics/DP4Pdf.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;
using namespace pybind11::literals;

void init_DP4Pdf(py::module &m) {
    py::class_<DPPdf, GooPdf>(m, "DPPdf")
        .def(py::init<std::string, std::vector<Observable>, DecayInfo4, GooPdf *>(),
             "n",
             "observables",
             "decay",
             "eff",
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>())
        .def(py::init<std::string, std::vector<Observable>, DecayInfo4, GooPdf *, unsigned int>(),
             "n",
             "observables",
             "decay",
             "eff",
             "MCeventsNorm",
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

                 py::array_t<fptype> pyparticles{{(size_t)4, 4 * numEvents}};
                 py::array_t<fptype> pyvariables{{(size_t)5, numEvents}};
                 py::array_t<fptype> pyweights{numEvents};
                 py::array_t<fptype> pyflags{numEvents};

                 for(int i = 0; i < 4; i++) {
                     for(int j = 0, k = 0; j < numEvents; j++, k = k + 4) {
                         pyparticles.mutable_at(i, k)     = (*(particles[i]))[j].get(0);
                         pyparticles.mutable_at(i, k + 1) = (*(particles[i]))[j].get(1);
                         pyparticles.mutable_at(i, k + 2) = (*(particles[i]))[j].get(2);
                         pyparticles.mutable_at(i, k + 3) = (*(particles[i]))[j].get(3);
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
}
