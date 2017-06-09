#include <pybind11/pybind11.h>

#include <goofit/PdfBase.h>
#include <goofit/FitManager.h>

namespace py = pybind11;
using namespace GooFit;

void init_FitManager(py::module &m) {
    py::class_<FitManager>(m, "FitManager")
        .def(py::init<PdfBase *>())
        // Can't directly wrap becase we (currently) don't want the return value in python
        .def("fit", [](FitManager &self) { self.fit(); });
}
