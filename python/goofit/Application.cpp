#include <goofit/Python.h>

#include <pybind11/iostream.h>

#include <goofit/Application.h>

#include <iostream>

using namespace GooFit;

void init_Application(py::module &m) {
    m.def("print_splash", &print_splash, "Print a splash screen", py::call_guard<py::scoped_ostream_redirect>());
    m.def("print_goofit_info",
          &print_goofit_info,
          "Print GooFit information",
          "gpuDevice"_a = 0,
          py::call_guard<py::scoped_ostream_redirect>())

        ;
}
