#include <goofit/Python.h>

#include <pybind11/stl.h>

#include <goofit/PDFs/physics/ThreeGaussResolutionSplice.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/physics/ThreeGaussResolutionSplice.h>

using namespace GooFit;

void init_ThreeGaussResolutionSplice(py::module &m) {
    py::class_<ThreeGaussResolutionSplice, MixingTimeResolution>(m, "ThreeGaussResolutionSplice")
        .def(py::init<Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      std::vector<Variable>,
                      std::vector<Variable>,
                      std::vector<Variable>,
                      std::vector<Variable>,
                      std::vector<Variable>>(),
             "cf"_a,
             "tf"_a,
             "cb"_a,
             "cs"_a,
             "tb"_a,
             "ts"_a,
             "ob"_a,
             "os"_a,
             "knots"_a,
             "a0"_a,
             "a1"_a,
             "a2"_a,
             "a3"_a)
        .def_static("help", []() { return HelpPrinter(ThreeGaussResolutionSplice_docs); });
}
