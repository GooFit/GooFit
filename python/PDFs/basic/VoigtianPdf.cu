#include <pybind11/pybind11.h>

#include <goofit/PDFs/basic/VoigtianPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_VoigtianPdf(py::module &m) {
    py::class_<VoigtianPdf, GooPdf>(m, "VoigtianPdf")
        .def(py::init<std::string, Variable *, Variable *, Variable *, Variable *>(),
             py::keep_alive<1,3>(),
             py::keep_alive<1,4>(),
             py::keep_alive<1,5>(),
             py::keep_alive<1,6>()
             );
}
