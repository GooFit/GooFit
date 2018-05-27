#include <goofit/Python.h>

#include <goofit/PDFs/basic/GaussianPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;

void init_GaussianPdf(py::module &m) {
    py::class_<GaussianPdf, GooPdf>(m, "GaussianPdf")
        .def(py::init<std::string, Observable, Variable, Variable>(), "name", "_x", "m", "s");
}
