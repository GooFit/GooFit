#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/physics/DP4Pdf.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/SpinFactors.h>

using namespace GooFit;
namespace py = pybind11;

void init_DP4Pdf(py::module &m) {
    py::class_<DPPdf, GooPdf>(m, "DPPdf")
        .def(py::init<std::string, std::vector<Variable *>, DecayInfo_DP *, GooPdf *>())
        .def(py::init<std::string, std::vector<Variable *>, DecayInfo_DP *, GooPdf *, unsigned int>())
        ;
}







