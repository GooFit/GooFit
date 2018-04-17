#include <pybind11/pybind11.h>

#include <goofit/PDFs/basic/BifurGaussPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_BifurGaussPdf(py::module &m) {
    py::class_<BifurGaussPdf, GooPdf>(m, "BifurGaussPdf")
        .def(py::init<std::string, Observable, Variable, Variable, Variable>(), "n", "_x", "m", "sL", "sR");
}
