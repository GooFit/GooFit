#include <pybind11/pybind11.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/basic/VoigtianPdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_VoigtianPdf(py::module &m) {
    py::class_<VoigtianPdf, GooPdf>(m, "VoigtianPdf")
        .def(py::init<std::string, Variable*, Variable*, Variable*, Variable*>())
        ;
}








