#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/PDFs/combine/EventWeightedAddPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_EventWeightedAddPdf(py::module &m) {
    py::class_<EventWeightedAddPdf, GooPdf>(m, "EventWeightedAddPdf")
        .def(py::init<std::string, std::vector<Observable>, std::vector<PdfBase *>>(), "n", "weights", "comps");
}
