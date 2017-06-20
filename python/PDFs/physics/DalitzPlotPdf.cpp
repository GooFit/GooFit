#include <pybind11/pybind11.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/physics/DalitzPlotPdf.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>

using namespace GooFit;
namespace py = pybind11;

void init_DalitzPlotPdf(py::module &m) {
    py::class_<DalitzPlotPdf, GooPdf>(m, "DalitzPlotPdf")
        .def(py::init<std::string, Variable *, Variable *, CountingVariable *, DecayInfo *, GooPdf *>())
        ;
}






