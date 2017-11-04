#include <pybind11/pybind11.h>

#include <goofit/PDFs/basic/NovosibirskPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_NovosibirskPdf(py::module &m) {
    py::class_<NovosibirskPdf, GooPdf>(m, "NovosibirskPdf")
        .def(py::init<std::string, Observable, Variable, Variable, Variable>(),
             py::keep_alive<1, 3>(),
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>(),
             py::keep_alive<1, 6>());
}
