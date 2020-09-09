#include <goofit/Python.h>

#include <pybind11/iostream.h>

#include <goofit/Application.h>

#include <iostream>

using namespace GooFit;

void init_Application(py::module &m) {
    m.def("print_splash", &print_splash, "Print a splash screen", py::call_guard<py::scoped_ostream_redirect>());

    m.def(
        "goofit_info",
        [](int gpuDev) { return goofit_info_version() + "\n" + goofit_info_device(gpuDev); },
        "Get GooFit information",
        "gpuDevice"_a = 0);

    m.def(
        "print_goofit_info",
        [](int gpuDev) { py::print(goofit_info_version() + "\n" + goofit_info_device(gpuDev)); },
        "Print GooFit information (same as goofit_info, kept for legacy code. Please use print(goofit_info) instead.)",
        "gpuDevice"_a = 0);

    m.def("set_floating_exceptions",
          &Application::set_floating_exceptions,
          "Enable floating point exceptions as Exceptions");

    m.def(
        "__main__", []() { print_goofit_info(0); }, py::call_guard<py::scoped_ostream_redirect>());
}
