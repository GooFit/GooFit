#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/FitControl.h>
#include <goofit/PDFs/GooPdf.h>

using namespace GooFit;
namespace py = pybind11;

using namespace pybind11::literals;

void init_GooPdf(py::module &m) {
    py::class_<GooPdf, PdfBase>(m, "GooPdf")
        .def("normalize", &GooPdf::normalize)
        .def("copyParams", (void (GooPdf::*)()) & GooPdf::copyParams)
        .def("normalize", &GooPdf::copyNormFactors)
        .def("makeGrid", &GooPdf::makeGrid)
        .def("evaluateAtPoints", &GooPdf::evaluateAtPoints)
        .def("setParameterConstantness", &GooPdf::setParameterConstantness)
        .def("evaluatePdf",
             [](GooPdf &self, Observable &var) {
                 auto grid     = self.makeGrid();
                 auto old_data = self.getData();
                 self.setData(&grid);
                 auto retval = self.evaluateAtPoints(var);
                 if(old_data != nullptr)
                     self.setData(old_data);
                 return py::make_tuple(grid, retval);
             },
             R"raw(
                Run makeGrid, set data, evaluateAtPoints, then recover original data.
                )raw")
        .def("setFitControl", &GooPdf::setFitControl, "Set a fit control.", "fit_control"_a)
        .def("getCompProbsAtDataPoints",
             &GooPdf::getCompProbsAtDataPoints,
             "Returns the probability at the current data points, supports multidimensional datasets.");
    ;
}
