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
        .def("setGenerationOffset", [](TDDP4& self, int off){})

        .def("GenerateSig",[](TDDP4 &self, int numEvents){
            mcbooster::ParticlesSet_h particles; // typedef for std::vector<Particles_h *>,
            mcbooster::VariableSet_h variables;
            mcbooster::RealVector_h weights;
            mcbooster::BoolVector_h flags;

            std::tie(particles, variables, weights, flags) = self.GenerateSig(numEvents);

            return std::make_tuple(particles, variables, weights, flags);
         }
    );
}


