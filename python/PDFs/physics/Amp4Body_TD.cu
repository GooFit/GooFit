#include <goofit/Python.h>

#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <goofit/DataSet.h>
#include <goofit/PDFs/physics/Amp4Body_TD.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/MixingTimeResolution.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

using namespace pybind11::literals;

void init_Amp4Body_TD(py::module &m) {
    py::class_<Amp4Body_TD, Amp4BodyBase> cls(m, "Amp4Body_TD");
    cls.def(py::init<std::string,
                     std::vector<Observable>,
                     DecayInfo4t,
                     MixingTimeResolution *,
                     GooPdf *,
                     Observable *,
                     long,
                     unsigned int>(),
            "n"_a,
            "observables"_a,
            "decay"_a,
            "r"_a,
            "eff"_a,
            "mistag"_a       = nullptr,
            "normSeed"_a     = 0,
            "MCeventsNorm"_a = 5e6,
            py::keep_alive<1, 4>(),
            py::keep_alive<1, 5>(),
            py::keep_alive<1, 6>(),
            py::keep_alive<1, 7>())

        .def("GenerateSig",
             [](Amp4Body_TD &self, size_t numEvents) {
                 mcbooster::ParticlesSet_h particles; // vector of pointers to vectors of 4R
                 mcbooster::VariableSet_h variables;  // vector of pointers to vectors of Grealt
                 mcbooster::RealVector_h weights;     // vector of greal t
                 mcbooster::BoolVector_h flags;       // vector of gboolt

                 std::tie(particles, variables, weights, flags) = self.GenerateSig(numEvents);

                 size_t nAcc = weights.size();

                 py::array_t<fptype> pyparticles{{(size_t)4, nAcc, (size_t)4}};
                 py::array_t<fptype> pyvariables{{(size_t)6, nAcc}};
                 py::array_t<fptype> pyweights{static_cast<py::ssize_t>(nAcc)};
                 py::array_t<bool> pyflags{static_cast<py::ssize_t>(nAcc)};

                 for(int i = 0; i < 4; i++) {
                     for(int j = 0; j < nAcc; j++) {
                         for(int k = 0; k < 4; k++) {
                             pyparticles.mutable_at(i, j, k) = (*(particles[i]))[j].get(k);
                         }
                     }
                 }

                 for(int i = 0; i < 6; i++) {
                     for(int j = 0; j < nAcc; j++) {
                         pyvariables.mutable_at(i, j) = (*(variables[i]))[j];
                     }
                 }

                 for(int i = 0; i < nAcc; i++) {
                     pyweights.mutable_at(i) = weights[i];
                 }

                 for(int i = 0; i < nAcc; i++) {
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
        .def("setGenerationOffset", &Amp4Body_TD::setGenerationOffset, "off"_a)

        ;

    m.attr("TDDP4") = cls;
}
