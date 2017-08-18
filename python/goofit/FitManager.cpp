#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>
#include <iostream>

#include <goofit/PdfBase.h>
#include <goofit/FitManager.h>

namespace py = pybind11;
using namespace GooFit;

void init_FitManager(py::module &m) {
    py::class_<FitManager>(m, "FitManager")
        .def(py::init<PdfBase *>())
        // Can't directly wrap becase we (currently) don't want the return value in python
        .def("fit", [](FitManager &self) {
                py::scoped_output_redirect redir(
                    std::cout, py::module::import("sys").attr("stdout")
                );
                self.fit();
            });
}
