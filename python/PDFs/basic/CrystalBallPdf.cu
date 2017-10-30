#include <pybind11/pybind11.h>

#include <goofit/PDFs/basic/CrystalBallPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_CrystalBallPdf(py::module &m) {
    py::class_<CrystalBallPdf, GooPdf>(m, "CrystalBallPdf")
        .def(py::init<std::string, Variable *, Variable *, Variable *, Variable *>(),
             py::keep_alive<1,3>(),
             py::keep_alive<1,4>(),
             py::keep_alive<1,5>(),
             py::keep_alive<1,6>())
        .def(py::init<std::string, Variable *, Variable *, Variable *, Variable *, Variable *>());
}
