
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/PDFs/GooPdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_GooPdf(py::module &m) {
    py::class_<GooPdf, PdfBase>(m, "GooPdf")
        .def("makeGrid", &GooPdf::makeGrid)
        .def("evaluateAtPoints", &GooPdf::evaluateAtPoints)
        .def("evaluatePdf", [](GooPdf &self, Variable &var){
                auto grid = self.makeGrid();
                auto old_data = self.getData();
                self.setData(&grid);
                auto retval = self.evaluateAtPoints(&var);
                if(old_data != nullptr)
                    self.setData(old_data);
                return py::make_tuple(grid, retval);
                },
                R"raw(
                Run makeGrid, set data, evaluateAtPoints, then recover original data.
                )raw")
    ;
}
