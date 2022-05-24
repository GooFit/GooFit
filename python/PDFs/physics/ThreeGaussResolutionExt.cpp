#include <goofit/Python.h>

#include <pybind11/stl.h>

#include <goofit/PDFs/physics/ThreeGaussResolutionExt.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/physics/ThreeGaussResolutionExt.h>

using namespace GooFit;

void init_ThreeGaussResolutionExt(py::module &m) {
    py::class_<ThreeGaussResolutionExt, MixingTimeResolution>(m, "ThreeGaussResolutionExt")
        .def(py::init<Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable>(),
             "cf"_a,
             "tf"_a,
             "cb"_a,
             "cs"_a,
             "tb"_a,
             "ts"_a,
             "ob"_a,
             "os"_a,
             "sb_low"_a,
             "sb_high"_a,
             "Tthres"_a,
             "constantC"_a)
        .def_static("help", []() { return HelpPrinter(ThreeGaussResolutionExt_docs); });
}
