#include <pybind11/pybind11.h>

#include <goofit/PDFs/basic/BifurGaussPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_BifurGaussPdf(py::module &m) {
    py::class_<BifurGaussPdf, GooPdf>(m, "BifurGaussPdf")
        .def(py::init<std::string, Observable, Variable, Variable, Variable>(),
             py::keep_alive<1, 3>(),
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>(),
             py::keep_alive<1, 6>());
}
