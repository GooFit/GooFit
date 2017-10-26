#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/DalitzPlotPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

using namespace pybind11::literals;

void init_DalitzPlotPdf(py::module &m) {
    py::class_<DalitzPlotPdf, GooPdf>(m, "DalitzPlotPdf")
        .def(py::init<std::string, Variable *, Variable *, CountingVariable *, DecayInfo *, GooPdf *>())
        .def("setDataSize", &DalitzPlotPdf::setDataSize, "dataSize"_a, "evtSize"_a = 3 )

        ;
}

