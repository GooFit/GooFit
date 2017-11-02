#include <pybind11/pybind11.h>

#include <goofit/PDFs/basic/CorrGaussianPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_CorrGaussianPdf(py::module &m) {
    py::class_<CorrGaussianPdf, GooPdf>(m, "CorrGaussianPdf")
        .def(py::init<std::string,
                      Observable,
                      Observable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable>(),
             py::keep_alive<1, 3>(),
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>(),
             py::keep_alive<1, 6>(),
             py::keep_alive<1, 7>(),
             py::keep_alive<1, 8>(),
             py::keep_alive<1, 9>());
}
