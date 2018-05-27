#include <goofit/Python.h>

#include <goofit/PDFs/basic/VoigtianPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;

void init_VoigtianPdf(py::module &m) {
    py::class_<VoigtianPdf, GooPdf>(m, "VoigtianPdf")
        .def(py::init<std::string, Observable, Variable, Variable, Variable>(), "n", "_x", "m", "s", "w");
}
