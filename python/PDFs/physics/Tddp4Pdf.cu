#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/physics/Tddp4Pdf.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/MixingTimeResolution_Aux.h>
#include <goofit/PDFs/physics/SpinFactors.h>

using namespace GooFit;
namespace py = pybind11;

void init_Tddp4Pdf(py::module &m) {
    py::class_<TDDP4, GooPdf>(m, "TDDP4")
        .def(py::init<std::string, std::vector<Variable *>, DecayInfo_DP *, MixingTimeResolution *, GooPdf*>())
        .def(py::init<std::string, std::vector<Variable *>, DecayInfo_DP *, MixingTimeResolution *, GooPdf*, Variable *>())
        .def(py::init<std::string, std::vector<Variable *>, DecayInfo_DP *, MixingTimeResolution *, GooPdf*, Variable *, unsigned int>())
        ;
}







