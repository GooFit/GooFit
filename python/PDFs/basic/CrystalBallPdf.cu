#include <pybind11/pybind11.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/basic/CrystalBallPdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_CrystalBallPdf(py::module &m) {
    py::class_<CrystalBallPdf, GooPdf>(m, "CrystalBallPdf")
        .def(py::init<std::string, Variable*, Variable*, Variable*, Variable*>())
        .def(py::init<std::string, Variable*, Variable*, Variable*, Variable*, Variable*>())
        ;
}


