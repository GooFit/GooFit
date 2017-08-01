#include <pybind11/pybind11.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/basic/ExpGausPdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_ExpGausPdf(py::module &m) {
    py::class_<ExpGausPdf, GooPdf>(m, "ExpGausPdf")
        .def(py::init<std::string, Variable*, Variable*, Variable*, Variable*>())
        ;
}


