#include <pybind11/pybind11.h>

#include <goofit/PDFs/basic/CorrGaussianPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_CorrGaussianPdf(py::module &m) {
    py::class_<CorrGaussianPdf, GooPdf>(m, "CorrGaussianPdf")
        .def(py::init<std::string, Observable, Observable, Variable, Variable, Variable, Variable, Variable>(),
             "n",
             "_x",
             "_y",
             "mean1",
             "sigma1",
             "mean2",
             "sigma2",
             "correlation");
}
