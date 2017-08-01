#include <pybind11/pybind11.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/physics/IncoherentSumPdf.h>
#include <goofit/PDFs/physics/TddpPdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_IncoherentSumPdf(py::module &m) {
    py::class_<IncoherentSumPdf, GooPdf>(m, "IncoherentSumPdf")
        .def(py::init<std::string, Variable *, Variable *, CountingVariable *, DecayInfo *, GooPdf *>())
        ;
}







