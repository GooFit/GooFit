#include <pybind11/pybind11.h>

#include <goofit/PDFs/basic/ArgusPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_ArgusPdf(py::module &m) {
    py::class_<ArgusPdf, GooPdf>(m, "ArgusPdf")
<<<<<<< HEAD
        .def(py::init<std::string, Variable *, Variable *, Variable *, bool>(),
             "n",
             "_x",
             "m",
             "s",
             "upper",
             py::keep_alive<1, 3>(),
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>())
        .def(py::init<std::string, Variable *, Variable *, Variable *, bool, Variable *>(),
             "n",
             "_x",
             "m",
             "s",
             "upper",
             "power",
             py::keep_alive<1, 3>(),
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>());
=======
        .def(py::init<std::string, Observable, Variable, Variable, bool>())
        .def(py::init<std::string, Observable, Variable, Variable, bool, Variable>());
>>>>>>> 3121d83c81449ab03f57f2444d77409085bdb7ea
}
