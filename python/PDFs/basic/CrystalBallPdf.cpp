#include <pybind11/pybind11.h>

#include <goofit/PDFs/basic/CrystalBallPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_CrystalBallPdf(py::module &m) {
    py::class_<CrystalBallPdf, GooPdf>(m, "CrystalBallPdf")
        .def(py::init<std::string, Observable, Variable, Variable, Variable>(), "n", "_x", "m", "s", "a")
        .def(py::init<std::string, Observable, Variable, Variable, Variable, Variable>(),
             "n",
             "_x",
             "m",
             "s",
             "a",
             "power");
}
