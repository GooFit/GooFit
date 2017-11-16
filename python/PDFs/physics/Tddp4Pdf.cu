#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/DataSet.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/MixingTimeResolution_Aux.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/PDFs/physics/Tddp4Pdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

using namespace pybind11::literals;

void init_Tddp4Pdf(py::module &m) {
    py::class_<TDDP4, GooPdf>(m, "TDDP4")
        .def(py::init<std::string,
                      std::vector<Observable>,
                      DecayInfo4t,
                      MixingTimeResolution *,
                      GooPdf *,
                      Observable *,
                      unsigned int>(),
             "n"_a,
             "observables"_a,
             "decay"_a,
             "r"_a,
             "eff"_a,
             "mistag"_a       = nullptr,
             "MCeventsNorm"_a = 5e6,
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>(),
             py::keep_alive<1, 6>(),
             py::keep_alive<1, 7>())

        .def("GenerateSig",
             [](TDDP4 &self, size_t numEvents) {
                 mcbooster::ParticlesSet_h particles; // vector of pointers to vectors of 4R
                 mcbooster::VariableSet_h variables;  // vector of pointers to vectors of Grealt
                 mcbooster::RealVector_h weights;     // vector of greal t
                 mcbooster::BoolVector_h flags;       // vector of gboolt

                 std::tie(particles, variables, weights, flags) = self.GenerateSig(numEvents);

                 py::array_t<fptype> pyparticles{{(size_t)4 * 4, numEvents}};
                 py::array_t<fptype> pyvariables{{(size_t)6, numEvents}};
                 py::array_t<fptype> pyweights{numEvents};
                 py::array_t<bool> pyflags{numEvents};

                 for(int i = 0; i < 4; i++) {
                     for(int j = 0; j < weights.size(); j++) {
                         pyparticles.mutable_at(i * 4, j)     = (*(particles[i]))[j].get(0);
                         pyparticles.mutable_at(i * 4 + 1, j) = (*(particles[i]))[j].get(1);
                         pyparticles.mutable_at(i * 4 + 2, j) = (*(particles[i]))[j].get(2);
                         pyparticles.mutable_at(i * 4 + 3, j) = (*(particles[i]))[j].get(3);
                     }
                 }

                 for(int i = 0; i < 6; i++) {
                     for(int j = 0; j < weights.size(); j++) {
                         pyvariables.mutable_at(i, j) = (*(variables[i]))[j];
                     }
                 }

                 for(int i = 0; i < weights.size(); i++) {
                     pyweights.mutable_at(i) = weights[i];
                 }

                 for(int i = 0; i < weights.size(); i++) {
                     pyflags.mutable_at(i) = flags[i];
                 }
                 delete variables[0];
                 delete variables[1];
                 delete variables[2];
                 delete variables[3];
                 delete variables[4];
                 delete variables[5];

                 delete particles[0];
                 delete particles[1];
                 delete particles[2];
                 delete particles[3];

                 return std::make_tuple(pyparticles, pyvariables, pyweights, pyflags);
             })
        .def("setGenerationOffset", &TDDP4::setGenerationOffset, "off"_a);
}
