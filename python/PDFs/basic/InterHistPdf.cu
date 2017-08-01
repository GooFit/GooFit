#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/basic/InterHistPdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_InterHistPdf(py::module &m) {
    py::class_<InterHistPdf, GooPdf>(m, "InterHistPdf")
        .def(py::init<std::string, BinnedDataSet *, std::vector<Variable *>, std::vector<Variable *>>())
        ;
}



