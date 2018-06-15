#include <goofit/Python.h>

#include <pybind11/stl.h>

#include <goofit/PDFs/physics/TruthResolution.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/physics/TruthResolution.h>

using namespace GooFit;

void init_TruthResolution(py::module &m) {
    py::class_<TruthResolution, MixingTimeResolution>(m, "TruthResolution")
        .def(py::init<>())
        .def_static("help", []() { return HelpPrinter(TruthResolution_docs); })

        ;
}
