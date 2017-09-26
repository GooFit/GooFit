#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/PDFs/physics/DP4Pdf.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;
using namespace pybind11::literals;

void init_DP4Pdf(py::module &m) {
    py::class_<DPPdf, GooPdf>(m, "DPPdf")
        .def(py::init<std::string,
                      std::vector<Variable *>,
                      DecayInfo_DP *,
                      GooPdf *>())
        .def(py::init<std::string,
                      std::vector<Variable *>,
                      DecayInfo_DP *,
                      GooPdf *,
                      unsigned int>())

        .def("normalize", &DPPdf::normalize)
        .def("setDataSize", &DPPdf::setDataSize, "dataSize"_a, "evtSize"_a = 6)
        .def("setForceIntegrals", &DPPdf::setForceIntegrals, "fset"_a = true)
        .def("getMCevents", &DPPdf::getMCevents)
        .def("setGenerationOffset", &DPPdf::setGenerationOffset, "off"_a)

        ;
}
