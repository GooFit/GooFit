#include <goofit/Python.h>

#include <goofit/PDFs/basic/NovosibirskPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;

void init_NovosibirskPdf(py::module &m) {
    py::class_<NovosibirskPdf, GooPdf>(m, "NovosibirskPdf")
        .def(py::init<std::string, Observable, Variable, Variable, Variable>(), "n", "_x", "m", "s", "t");
}
