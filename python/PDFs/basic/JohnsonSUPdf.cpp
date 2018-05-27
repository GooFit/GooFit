#include <goofit/Python.h>

#include <goofit/PDFs/basic/JohnsonSUPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;

void init_JohnsonSUPdf(py::module &m) {
    py::class_<JohnsonSUPdf, GooPdf>(m, "JohnsonSUPdf")
        .def(
            py::init<std::string, Observable, Variable, Variable, Variable, Variable>(), "n", "_x", "m", "s", "g", "d");
}
