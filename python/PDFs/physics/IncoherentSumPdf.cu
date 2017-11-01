#include <pybind11/pybind11.h>

#include <goofit/PDFs/physics/IncoherentSumPdf.h>
#include <goofit/PDFs/physics/TddpPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_IncoherentSumPdf(py::module &m) {
    py::class_<IncoherentSumPdf, GooPdf>(m, "IncoherentSumPdf")
        .def(py::init<std::string, Observable *, Observable *, EventNumber *, DecayInfo *, GooPdf *>(),
             py::keep_alive<1, 3>(),
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>(),
             py::keep_alive<1, 6>(),
             py::keep_alive<1, 7>());
}
