#include <pybind11/pybind11.h>

#include <goofit/PDFs/basic/TrigThresholdPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;
using namespace pybind11::literals;

void init_TrigThresholdPdf(py::module &m) {
    py::class_<TrigThresholdPdf, GooPdf>(m, "TrigThresholdPdf")
        .def(py::init<std::string, Observable *, Variable *, Variable *, Variable *, bool>(),
             "n"_a, "x"_a, "thresh"_a, "trigConst"_a, "linConst"_a, "upper"_a = true,
             py::keep_alive<1, 3>(),
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>(),
             py::keep_alive<1, 6>())
        .def(py::init<std::string, Observable *, Observable *, Variable *, Variable *, Variable *, Variable *, bool>(),
             "n"_a, "x"_a, "y"_a,"thresh"_a, "trigConst"_a, "linConst"_a, "massConstant"_a, "upper"_a,
             py::keep_alive<1, 3>(),
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>(),
             py::keep_alive<1, 6>(),
             py::keep_alive<1, 7>(),
             py::keep_alive<1, 8>());
}
