#include <pybind11/pybind11.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/basic/CorrGaussianPdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_CorrGaussianPdf(py::module &m) {
    py::class_<CorrGaussianPdf, GooPdf>(m, "CorrGaussianPdf")
        .def(py::init<std::string, Variable*, Variable*, Variable*, Variable*, Variable*, Variable*, Variable*>())
        ;
}


