#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/PDFs/combine/MappedPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_MappedPdf(py::module &m) {
    py::class_<MappedPdf, GooPdf>(m, "MappedPdf")
        .def(py::init<std::string, GooPdf *, std::vector<GooPdf *> &>(),
             "n",
             "m",
             "t",
             py::keep_alive<1, 3>(),
             py::keep_alive<1, 4>());
}
