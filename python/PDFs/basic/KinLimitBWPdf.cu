#include <pybind11/pybind11.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/basic/KinLimitBWPdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_KinLimitBWPdf(py::module &m) {
    py::class_<KinLimitBWPdf, GooPdf>(m, "KinLimitBWPdf")
        .def(py::init<std::string, Variable*, Variable*, Variable*>())
        ;
}




