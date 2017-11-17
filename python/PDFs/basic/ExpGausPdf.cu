#include <pybind11/pybind11.h>

#include <goofit/PDFs/basic/ExpGausPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_ExpGausPdf(py::module &m) {
    py::class_<ExpGausPdf, GooPdf>(m, "ExpGausPdf")
        .def(py::init<std::string, Observable, Variable, Variable, Variable>(), "n", "_x", "m", "s", "t");
}
