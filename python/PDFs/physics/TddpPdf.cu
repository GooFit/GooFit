#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/physics/TddpPdf.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/MixingTimeResolution_Aux.h>

using namespace GooFit;
namespace py = pybind11;

void init_TddpPdf(py::module &m) {
    py::class_<TddpPdf, GooPdf>(m, "TddpPdf")
        .def(py::init<std::string, Variable *, Variable *,Variable *, Variable *, CountingVariable *, DecayInfo *, MixingTimeResolution *, GooPdf *>())
        .def(py::init<std::string, Variable *, Variable *,Variable *, Variable *, CountingVariable *, DecayInfo *, MixingTimeResolution *, GooPdf *, Variable *>())
        .def(py::init<std::string, Variable *, Variable *,Variable *, Variable *, CountingVariable *, DecayInfo *, std::vector<MixingTimeResolution *>, GooPdf *, Variable *>())
        .def(py::init<std::string, Variable *, Variable *,Variable *, Variable *, CountingVariable *, DecayInfo *, std::vector<MixingTimeResolution *>, GooPdf *, Variable *, Variable *>())
        ;
}








