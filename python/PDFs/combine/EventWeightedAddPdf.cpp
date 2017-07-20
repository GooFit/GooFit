#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/combine/EventWeightedAddPdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_EventWeightedAddPdf(py::module &m) {
    py::class_<EventWeightedAddPdf, GooPdf>(m, "EventWeightedAddPdf")
        .def(py::init<std::string, std::vector<Variable *>, std::vector<PdfBase *>>())
        ;
}





