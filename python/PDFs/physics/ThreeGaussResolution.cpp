#include <goofit/Python.h>

#include <pybind11/stl.h>

#include <goofit/PDFs/physics/ThreeGaussResolution.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/physics/ThreeGaussResolution.h>

using namespace GooFit;

void init_ThreeGaussResolution(py::module &m) {
    py::class_<ThreeGaussResolution, MixingTimeResolution>(m, "ThreeGaussResolution")
        .def(py::init<Variable, Variable, Variable, Variable, Variable, Variable, Variable, Variable, Variable>(),
             "cf"_a,
             "tf"_a,
             "cb"_a,
             "cs"_a,
             "tb"_a,
             "ts"_a,
             "ob"_a,
             "os"_a,
             "sb"_a)

        .def_static("help", []() { return HelpPrinter(ThreeGaussResolution_docs); });
}
