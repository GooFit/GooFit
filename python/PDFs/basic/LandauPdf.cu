#include <pybind11/pybind11.h>

#include <goofit/PDFs/basic/LandauPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_LandauPdf(py::module &m) {
    py::class_<LandauPdf, GooPdf>(m, "LandauPdf")
        .def(py::init<std::string, Observable, Variable, Variable>(), "n", "_x", "mpv", "sigma");
}
