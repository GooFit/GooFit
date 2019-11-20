#include <goofit/Python.h>

#include <pybind11/complex.h>
#include <pybind11/stl.h>

#include <goofit/PDFs/physics/SquareDalitzEffPdf.h>

using namespace GooFit;

void init_SquareDalitzEffPdf(py::module &m) {
    py::class_<SquareDalitzEffPdf, GooPdf>(m, "SquareDalitzEffPdf")
        .def(py::init<std::string,
                      Observable,
                      Observable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable>(),
             "name"_a,
             "m12"_a,
             "m13"_a,
             "c0"_a,
             "c1"_a,
             "c2"_a,
             "c3"_a,
             "c4"_a,
             "c5"_a,
             "c6"_a,
             "c7"_a);
}
