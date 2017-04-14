#include <pybind11/pybind11.h>
#include <goofit/PdfBase.h>
#include <goofit/FitManager.h>

namespace py = pybind11;

void init_FitManager(py::module &m) {
    py::class_<FitManager>(m, "FitManager")
        .def(py::init<PdfBase*>())
        .def("fit", &FitManager::fit)
        ;


}
