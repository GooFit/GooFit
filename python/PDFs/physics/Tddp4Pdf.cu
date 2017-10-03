#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

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
                      std::vector<Variable *>,
                      DecayInfo_DP *,
                      MixingTimeResolution *,
                      GooPdf *>())
        .def(py::init<std::string,
                      std::vector<Variable *>,
                      DecayInfo_DP *,
                      MixingTimeResolution *,
                      GooPdf *,
                      Variable *>())     
        .def(py::init<std::string,
                      std::vector<Variable *>,
                      DecayInfo_DP *,
                      MixingTimeResolution *,
                      GooPdf *,
                      Variable *,
                      unsigned int>())
        .def("setGenerationOffset", &TDDP4::setGenerationOffset, "off"_a)
        /*
        .def("GenerateSig",[](TDDP4 &self, int numEvents){
            mcbooster::ParticlesSet_h particles; // typedef for std::vector<Particles_h *>,
            mcbooster::VariableSet_h variables;
            mcbooster::RealVector_h weights;
            mcbooster::BoolVector_h flags;

            std::tie(particles, variables, weights, flags) = self.GenerateSig(numEvents);

            py::array_t<fptype> pyparticles{{16, numEvents}};

            for(int i = 0; i < particles.size()*4; i++)
                for(int j = 0; j < numEvents; j++)
                // (*(particles[i/4]))[j].get(i%4)
                    result.mutable_at(i, j) = instance.getValue(instance.getVariables().at(i), j);

            return std::make_tuple(pyparticles, pyvariables, pyweights, pyflags);
         })
         */
    ;
}


