#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/PDFs/combine/AddPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_AddPdf(py::module &m) {
    py::class_<AddPdf, GooPdf>(m, "AddPdf")
        .def(py::init<std::string, std::vector<Variable>, std::vector<PdfBase *>>(), "n", "weights", "comps")
        .def(py::init<std::string, Variable, PdfBase *, PdfBase *>(),
             "n",
             "frac1",
             "func1",
             "func2",
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>());
}
