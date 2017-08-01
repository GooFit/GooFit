#include <pybind11/pybind11.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/basic/TrigThresholdPdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_TrigThresholdPdf(py::module &m) {
    py::class_<TrigThresholdPdf, GooPdf>(m, "TrigThresholdPdf")
        .def(py::init<std::string, Variable*, Variable*, Variable*, Variable*>())
        .def(py::init<std::string, Variable*, Variable*, Variable*, Variable*,  bool>())
        .def(py::init<std::string, Variable*, Variable*, Variable*, Variable*, Variable*, Variable*, bool>())
        ;
}







