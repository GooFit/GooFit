#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/MixingTimeResolution_Aux.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/PDFs/physics/Tddp4Pdf.h>
#include <goofit/Variable.h>
#include <goofit/DataSet.h>

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

        .def("GenerateSig",[](TDDP4 &self, int numEvents){
            mcbooster::ParticlesSet_h particles; // typedef for std::vector<Particles_h *>,
            mcbooster::VariableSet_h variables;
            mcbooster::RealVector_h weights;
            mcbooster::BoolVector_h flags;

            std::tie(particles, variables, weights, flags) = self.GenerateSig(numEvents);

            py::array_t<fptype> pyparticles{{4, 4*numEvents}};
            py::array_t<fptype> pyvariables{{6, numEvents}};
            py::array_t<fptype> pyweights{numEvents};
            py::array_t<fptype> pyflags{numEvents};

            for(int i = 0; i < 4; i++){
                for(int j = 0, k = 0; j < numEvents; j++, k=k+4){
                    pyparticles.mutable_at(i, k) = (*(particles[i]))[j].get(0);
                    pyparticles.mutable_at(i, k+1) = (*(particles[i]))[j].get(1);
                    pyparticles.mutable_at(i, k+2) = (*(particles[i]))[j].get(2);
                    pyparticles.mutable_at(i, k+3) = (*(particles[i]))[j].get(3);
                }
            }

            for(int i = 0; i < 6; i++){
                for(int j=0; j < numEvents; j++){
                    pyvariables.mutable_at(i, j) = (*(variables[i]))[j];
                }
            }

            for(int i = 0; i < numEvents; i++){
                pyweights.mutable_at(i) = weights[i];
            }

            for(int i = 0; i < numEvents; i++){
                pyflags.mutable_at(i) = flags[i];
            }

            return std::make_tuple(pyparticles, pyvariables, pyweights, pyflags);
         })

    ;
}


